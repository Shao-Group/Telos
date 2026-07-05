"""
Stage I **BAM** feature extraction at candidate TSS/TES sites (:class:`~telos.candidates.extract.CandidateSite`).

Computes read density, coverage before/after the site, splice distances, soft-clip bias, gradient
sharpness, and k-mer statistics on clipped sequence. Supports sequential or multiprocess execution
with one pysam :class:`~pysam.AlignmentFile` per worker.
"""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import Pool, cpu_count
from pathlib import Path
from typing import Iterable
import csv
import math
import re
from collections import Counter, OrderedDict

import numpy as np
import pysam  # type: ignore

from telos.candidates.extract import CandidateSite

# Open once per worker process (set by Pool initializer).
_STAGE1_POOL_BAM: object | None = None


def _stage1_pool_init(bam_path: str) -> None:
    """Multiprocessing initializer: open ``bam_path`` once per worker and store in ``_STAGE1_POOL_BAM``."""
    global _STAGE1_POOL_BAM
    _STAGE1_POOL_BAM = pysam.AlignmentFile(bam_path, "rb")


@dataclass(frozen=True)
class Stage1FeatureConfig:
    """
    Numeric hyperparameters for windows, MAPQ cutoff, and k-mer list.

    Populated from YAML ``stage1.feature_extraction`` and ``stage1.features`` via
    :func:`telos.pipeline_core.build_stage1_runtime_config`.
    """

    window_size: int = 100
    density_window: int = 100
    coverage_window: int = 100
    soft_clip_window: int = 10
    min_mapq: int = 10
    splice_site_window: int = 300
    gradient_analysis_range: int = 100
    extended_window: int = 1000
    k_values: tuple[int, ...] = (3,)


def compute_strand(read: pysam.AlignedSegment) -> str:
    """
    Infer read strand using ``ts`` (transcription strand) BAM tag when present, else ``is_reverse``.

    If ``ts`` is ``+`` or ``-``, flips relative to read orientation so the returned value is genomic strand.
    Returns ``.`` only when ``ts`` exists but is not ``+``/``-``.
    """
    try:
        ts = read.get_tag("ts")
        if ts not in ["+", "-"]:
            return "."
        if read.is_reverse:
            return "-" if ts == "+" else "+"
        return "+" if ts == "+" else "-"
    except KeyError:
        return "-" if read.is_reverse else "+"


def calculate_entropy(positions: np.ndarray) -> float:
    """Shannon entropy (base 2) of discrete positions treated as categorical draws; ``0`` if empty."""
    if len(positions) == 0:
        return 0.0
    count = Counter(positions.tolist())
    total = sum(count.values())
    return float(-sum((freq / total) * math.log2(freq / total) for freq in count.values()))


def read_start_end_entropy(
    start_positions: np.ndarray, end_positions: np.ndarray
) -> tuple[float, float]:
    """Return ``(entropy(starts), entropy(ends))``"""
    return calculate_entropy(start_positions), calculate_entropy(end_positions)


def _read_overlaps_half_open(rs: int, re: int, start: int, end: int) -> bool:
    """True if alignment [rs, re) overlaps BAM fetch interval [start, end)."""
    return re > start and rs < end


def _accumulate_coverage_before_after(
    read: pysam.AlignedSegment, pos: int, cw: int, coverage_before: float, coverage_after: float
) -> tuple[float, float]:
    """
    Add overlap length (in bp) between read and ``[pos-cw, pos)`` (before) and ``[pos, pos+cw)`` (after).

    Uses half-open intervals on reference coordinates; requires mapped read with valid reference span.
    """
    reference_start = read.reference_start
    reference_end = read.reference_end
    if reference_start is None or reference_end is None:
        return coverage_before, coverage_after
    if pos - cw <= reference_start < pos or pos - cw < reference_end <= pos:
        overlap_start = max(reference_start, pos - cw)
        overlap_end = min(reference_end, pos)
        if overlap_end > overlap_start:
            coverage_before += overlap_end - overlap_start
    if pos <= reference_start < pos + cw or pos < reference_end <= pos + cw:
        overlap_start = max(reference_start, pos)
        overlap_end = min(reference_end, pos + cw)
        if overlap_end > overlap_start:
            coverage_after += overlap_end - overlap_start
    return coverage_before, coverage_after


def collect_narrow_reads_and_coverage(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos: int,
    cfg: Stage1FeatureConfig,
    strand: str,
) -> tuple[list[pysam.AlignedSegment], float, float, float]:
    """
    Fetch reads on ``chrom`` in a wide window, accumulate coverage before/after ``pos``, collect narrow window reads.

    Returns ``(narrow_reads, cov_before, cov_after, delta)`` where ``delta`` is oriented so ``+`` strand
    uses ``after - before`` and ``-`` strand flips sign so upstream/downstream semantics match transcription.
    """
    half = max(cfg.extended_window, cfg.window_size)
    fetch_start = max(0, pos - half)
    fetch_end = pos + half
    narrow_start = max(0, pos - cfg.window_size)
    narrow_end = pos + cfg.window_size

    coverage_before = 0.0
    coverage_after = 0.0
    narrow_reads: list[pysam.AlignedSegment] = []
    cw = cfg.coverage_window
    minq = cfg.min_mapq

    for read in bam.fetch(chrom, fetch_start, fetch_end):
        rs, re = read.reference_start, read.reference_end
        if rs is None or re is None:
            continue
        if not read.is_unmapped and read.mapping_quality >= minq:
            coverage_before, coverage_after = _accumulate_coverage_before_after(
                read, pos, cw, coverage_before, coverage_after
            )
        if _read_overlaps_half_open(rs, re, narrow_start, narrow_end):
            narrow_reads.append(read)

    delta_coverage = coverage_after - coverage_before
    if strand == "+":
        return narrow_reads, coverage_before, coverage_after, delta_coverage
    return narrow_reads, coverage_after, coverage_before, -delta_coverage


def nearest_splice_and_softclip_bias(
    reads: list[pysam.AlignedSegment], pos: int, cfg: Stage1FeatureConfig, strand: str
) -> tuple[int, float]:
    """
    From overlapping reads: minimum distance from ``pos`` to any N-op (splice) in MAPQ-passing alignments,
    and fraction of strand-matching reads with soft-clip within 5bp of ``pos`` (5' or 3' end by cigar).
    """
    min_distance = cfg.splice_site_window
    soft_clip_count = 0
    total_reads = 0
    minq = cfg.min_mapq

    for read in reads:
        if not read.is_unmapped and read.mapping_quality >= minq and read.cigartuples:
            curr_pos = read.reference_start
            for op, length in read.cigartuples:
                if op in (0, 7, 8):
                    curr_pos += length
                elif op == 3:
                    splice_start = curr_pos
                    splice_end = curr_pos + length
                    distance = min(abs(pos - splice_start), abs(pos - splice_end))
                    min_distance = min(min_distance, distance)
                    curr_pos += length
                else:
                    curr_pos += length

        read_strand = compute_strand(read)
        if (
            read.is_unmapped
            or read.mapping_quality < minq
            or not read.cigartuples
            or read_strand != strand
        ):
            continue
        total_reads += 1
        if read.cigartuples[0][0] == 4 and abs(read.reference_start - pos) <= 5:
            soft_clip_count += 1
        elif read.cigartuples[-1][0] == 4 and abs(read.reference_end - pos) <= 5:
            soft_clip_count += 1

    return int(min_distance), (soft_clip_count / total_reads if total_reads else 0.0)


def calculate_coverage_gradient_sharpness(
    cached_reads: list[pysam.AlignedSegment], pos: int, cfg: Stage1FeatureConfig, strand: str
) -> tuple[float, float, float]:
    """
    Bin coverage along a grid around ``pos``, take first differences, and summarize sharpness vs local mean.

    Returns ``(gradient_sharpness, max_abs_gradient, local_avg_coverage)``. For ``-`` strand, uses
    negated gradients for directional sharpness.
    """
    window_size = 10
    analysis_range = cfg.gradient_analysis_range
    positions = range(pos - analysis_range, pos + analysis_range + 1, window_size)
    coverages = []
    for window_start in positions:
        window_end = window_start + window_size
        coverage = 0
        for read in cached_reads:
            if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
                continue
            overlap_start = max(read.reference_start, window_start)
            overlap_end = min(read.reference_end, window_end)
            if overlap_end > overlap_start:
                coverage += overlap_end - overlap_start
        coverages.append(coverage / window_size)
    if len(coverages) < 3:
        return 0.0, 0.0, 0.0
    cov = np.array(coverages)
    gradients = np.diff(cov)
    max_gradient = float(np.max(np.abs(gradients))) if len(gradients) > 0 else 0.0
    local_avg_coverage = float(np.mean(cov)) if len(cov) > 0 else 1.0
    gradient_sharpness = max_gradient / max(local_avg_coverage, 1.0)
    if strand == "-":
        signed_gradients = -gradients
        max_signed_gradient = float(np.max(signed_gradients)) if len(signed_gradients) > 0 else 0.0
        if max_signed_gradient > 0:
            gradient_sharpness = max_signed_gradient / max(local_avg_coverage, 1.0)
    return gradient_sharpness, max_gradient, local_avg_coverage


def calculate_skewness(values: np.ndarray) -> float:
    """Third standardized moment; ``0`` if fewer than 3 values or zero std."""
    if len(values) < 3:
        return 0.0
    mean_val = float(np.mean(values))
    std_val = float(np.std(values))
    if std_val == 0:
        return 0.0
    skew = float(np.mean(((values - mean_val) ** 3))) / (std_val ** 3)
    return skew


def calculate_degradation_score(read_starts: np.ndarray, pos: int) -> float:
    """
    Negative slope of read-start counts across small windows around ``pos`` (higher = sharper drop-off).

    Bins ``[-50, 50]`` in steps of 10bp, fits least-squares slope of counts vs bin index, returns ``-slope``.
    """
    if len(read_starts) == 0:
        return 0.0
    window_size = 10
    windows = []
    for i in range(-50, 51, window_size):
        window_start = pos + i
        window_end = pos + i + window_size
        count = np.sum((read_starts >= window_start) & (read_starts < window_end))
        windows.append(int(count))
    if not windows or sum(windows) == 0:
        return 0.0
    x = list(range(len(windows)))
    y = windows
    n = len(x)
    slope = (n * sum(x[i] * y[i] for i in range(n)) - sum(x) * sum(y)) / max(
        n * sum(x[i] ** 2 for i in range(n)) - sum(x) ** 2, 1
    )
    return float(-slope)


def extract_stage1_stat_features(
    read_starts: np.ndarray,
    start_soft_clips: np.ndarray,
    end_soft_clips: np.ndarray,
    coverage_before: float,
    coverage_after: float,
    total_reads: int,
    pos: int,
) -> dict[str, float]:
    """
    Aggregate scalar features from read start positions, soft-clip arrays, and coverage totals.

    Includes coverage ratios, read-start variance/clustering, soft-clip sparsity and length moments,
    and normalization by coverage. Keys match Stage I training column names expected downstream.
    """
    features: dict[str, float] = {}
    features["coverage_ratio"] = coverage_after / max(coverage_before, 1.0)
    features["coverage_log_ratio"] = (
        math.log2(max(coverage_after, 1.0) / max(coverage_before, 1.0)) if coverage_before > 0 else 0.0
    )
    features["normalized_delta_coverage"] = (coverage_after - coverage_before) / max(
        coverage_after + coverage_before, 1.0
    )

    if len(read_starts) > 0:
        read_start_variance = float(np.var(read_starts))
        read_start_mean = float(np.mean(read_starts))
        features["read_start_variance"] = read_start_variance
        features["read_start_coefficient_variation"] = read_start_variance / max(read_start_mean, 1.0)
        distances = np.abs(read_starts - pos)
        features["weighted_read_start_density"] = float(np.sum(1 / np.maximum(distances, 1)))
        close_reads = int(np.sum(distances <= 5))
        features["read_start_clustering"] = close_reads / max(len(read_starts), 1)
        upstream_starts = int(np.sum(read_starts < pos))
        downstream_starts = int(np.sum(read_starts > pos))
        features["up_down_stream_ratio"] = upstream_starts / max(downstream_starts, 1)
        features["five_prime_degradation_score"] = calculate_degradation_score(read_starts, pos)
    else:
        features["read_start_variance"] = 0.0
        features["read_start_coefficient_variation"] = 0.0
        features["weighted_read_start_density"] = 0.0
        features["read_start_clustering"] = 0.0
        features["upstream_downstream_ratio"] = 0.0
        features["five_prime_degradation_score"] = 0.0

    total_soft_clips = len(start_soft_clips) + len(end_soft_clips)
    features["softclip_sparsity"] = total_soft_clips / max(total_reads, 1)
    if total_soft_clips > 0:
        all_soft = np.concatenate([start_soft_clips, end_soft_clips])
        features["softclip_length_variance"] = float(np.var(all_soft))
        features["softclip_length_skewness"] = float(calculate_skewness(all_soft))
    else:
        features["softclip_length_variance"] = 0.0
        features["softclip_length_skewness"] = 0.0

    features["norm_read_density"] = len(read_starts) / max(coverage_after, 1.0)
    features["coverage_efficiency"] = total_reads / max(coverage_before + coverage_after, 1.0)
    return features


def softclip_kmer_features(
    clipped_sequences: list[str], prefix: str = "", klist: tuple[int, ...] = (3,)
) -> dict[str, float]:
    """
    Composition and k-mer statistics over soft-clipped sequences near the site.

    Adds GC/purine content, max homopolymer runs per base, and :func:`_extract_optimized_kmer_features`
    for each ``k`` in ``klist``. Prefix keys when embedding in larger feature dicts.
    """
    if not clipped_sequences:
        clipped_sequences = [""]
    combined: dict[str, float] = {}
    processed = [seq.upper() for seq in clipped_sequences if seq]
    all_bases = "".join(processed)
    total_bases = len(all_bases)
    if total_bases > 0:
        a_count = all_bases.count("A")
        t_count = all_bases.count("T")
        g_count = all_bases.count("G")
        c_count = all_bases.count("C")
        combined[f"{prefix}gc_content"] = (g_count + c_count) / total_bases
        combined[f"{prefix}purine_ratio"] = (a_count + g_count) / total_bases
    else:
        combined[f"{prefix}gc_content"] = 0.0
        combined[f"{prefix}purine_ratio"] = 0.0

    for base in "ACGT":
        max_run = 0
        pattern = re.compile(f"{base}+")
        for seq in processed:
            matches = pattern.findall(seq)
            if matches:
                seq_max = max(len(m) for m in matches)
                max_run = max(max_run, seq_max)
        combined[f"{prefix}max_poly{base}"] = float(max_run)

    combined.update(_extract_optimized_kmer_features(processed, prefix, klist))
    return combined


def _extract_optimized_kmer_features(
    processed_sequences: list[str], prefix: str, klist: tuple[int, ...]
) -> dict[str, float]:
    """
    Count k-mers (k in ``klist``), fixed motif hits for k=3, diversity/entropy, and composition buckets.

    Only k-mers over ``ATGC`` are counted. Emits many ``{prefix}k{k}_*`` feature keys; zeros when no k-mers.
    """
    features: dict[str, float] = {}
    fixed_kmers_3 = {
        "ATG",
        "TAA",
        "TGA",
        "TAG",
        "GCC",
        "CGC",
        "GCG",
        "AAA",
        "TTT",
        "AAT",
        "GTG",
        "CAG",
        "GAG",
        "CCC",
        "GGG",
        "CCG",
    }
    kmer_data = {
        k: {
            "counts": Counter(),
            "fixed_counts": Counter(),
            "composition_counts": {
                "at_rich_kmers": 0,
                "gc_rich_kmers": 0,
                "purine_rich_kmers": 0,
                "repeat_kmers": 0,
            },
            "total_kmers": 0,
            "gc_kmers": 0,
        }
        for k in klist
    }
    for seq in processed_sequences:
        seq_len = len(seq)
        for k in klist:
            if seq_len < k:
                continue
            data = kmer_data[k]
            for i in range(seq_len - k + 1):
                kmer = seq[i : i + k]
                if all(base in "ATGC" for base in kmer):
                    data["counts"][kmer] += 1
                    data["total_kmers"] += 1
                    if k == 3 and kmer in fixed_kmers_3:
                        data["fixed_counts"][kmer] += 1
                    gc_count = kmer.count("G") + kmer.count("C")
                    a_count = kmer.count("A")
                    g_count = kmer.count("G")
                    gc_content = gc_count / k
                    if gc_content >= 0.5:
                        data["gc_kmers"] += 1
                    purine_content = (a_count + g_count) / k
                    if len(set(kmer)) == 1:
                        data["composition_counts"]["repeat_kmers"] += 1
                    elif gc_content > 0.6:
                        data["composition_counts"]["gc_rich_kmers"] += 1
                    elif gc_content < 0.4:
                        data["composition_counts"]["at_rich_kmers"] += 1
                    else:
                        pass
                    if purine_content > 0.6:
                        data["composition_counts"]["purine_rich_kmers"] += 1

    for k in klist:
        data = kmer_data[k]
        if k == 3:
            for kmer in fixed_kmers_3:
                features[f"{prefix}kmer_{kmer}"] = float(data["fixed_counts"].get(kmer, 0))
        total_kmers = data["total_kmers"]
        if total_kmers > 0:
            counts = data["counts"]
            unique_kmers = len(counts)
            diversity = unique_kmers / total_kmers
            most_frequent_count = max(counts.values()) if counts else 0
            most_frequent_ratio = most_frequent_count / total_kmers
            counts_array = np.array(list(counts.values()))
            probs = counts_array / total_kmers
            entropy = float(-np.sum(probs * np.log2(probs)))
            features[f"{prefix}k{k}_gc_kmers_ratio"] = data["gc_kmers"] / total_kmers
            features[f"{prefix}k{k}_kmer_diversity"] = diversity
            features[f"{prefix}k{k}_kmer_entropy"] = entropy
            features[f"{prefix}k{k}_most_frequent_kmer_ratio"] = most_frequent_ratio
            features[f"{prefix}k{k}_unique_kmers"] = float(unique_kmers)
            for comp_type, count in data["composition_counts"].items():
                features[f"{prefix}k{k}_{comp_type}"] = float(count)
        else:
            features[f"{prefix}k{k}_gc_kmers_ratio"] = 0.0
            features[f"{prefix}k{k}_kmer_diversity"] = 0.0
            features[f"{prefix}k{k}_kmer_entropy"] = 0.0
            features[f"{prefix}k{k}_most_frequent_kmer_ratio"] = 0.0
            features[f"{prefix}k{k}_unique_kmers"] = 0.0
            for comp_type in data["composition_counts"].keys():
                features[f"{prefix}k{k}_{comp_type}"] = 0.0
    return features


def compute_stage1_features_for_site(
    bam: pysam.AlignmentFile, site: CandidateSite, cfg: Stage1FeatureConfig
) -> dict[str, object]:
    """
    Full feature dict for one ``site``: fetch reads, compute scalars, append k-mer features, stable key order.

    Returns a plain dict with identity columns first, then sorted remaining keys (for reproducible TSV export).
    """
    chrom = site.chrom
    pos = int(site.position)
    strand = site.strand
    narrow_reads, coverage_before, coverage_after, delta_coverage = collect_narrow_reads_and_coverage(
        bam, chrom, pos, cfg, strand
    )

    read_data: dict[str, list] = {
        "starts": [],
        "ends": [],
        "start_soft_clips": [],
        "end_soft_clips": [],
        "map_quals": [],
        "strands": [],
        "start_clip_seqs": [],
        "end_clip_seqs": [],
    }
    valid_reads = []
    for read in narrow_reads:
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
            continue
        valid_reads.append(read)
        read_strand = compute_strand(read)
        left_soft_clip_length = 0
        right_soft_clip_length = 0
        left_soft_clip_seq = ""
        right_soft_clip_seq = ""
        if read.cigartuples and read.query_sequence:
            if read.cigartuples[0][0] == 4:
                left_soft_clip_length = read.cigartuples[0][1]
                left_soft_clip_seq = read.query_sequence[:left_soft_clip_length]
            if read.cigartuples[-1][0] == 4:
                right_soft_clip_length = read.cigartuples[-1][1]
                right_soft_clip_seq = read.query_sequence[-right_soft_clip_length:]
        if read_strand == "+":
            read_data["starts"].append(read.reference_start)
            read_data["ends"].append(read.reference_end)
            read_data["start_soft_clips"].append(left_soft_clip_length)
            read_data["end_soft_clips"].append(right_soft_clip_length)
            if left_soft_clip_seq and abs(read.reference_start - pos) <= cfg.soft_clip_window:
                read_data["start_clip_seqs"].append(left_soft_clip_seq)
            if right_soft_clip_seq and abs(read.reference_end - pos) <= cfg.soft_clip_window:
                read_data["end_clip_seqs"].append(right_soft_clip_seq)
        else:
            read_data["starts"].append(read.reference_end)
            read_data["ends"].append(read.reference_start)
            read_data["start_soft_clips"].append(right_soft_clip_length)
            read_data["end_soft_clips"].append(left_soft_clip_length)
            if right_soft_clip_seq and abs(read.reference_end - pos) <= cfg.soft_clip_window:
                read_data["start_clip_seqs"].append(right_soft_clip_seq)
            if left_soft_clip_seq and abs(read.reference_start - pos) <= cfg.soft_clip_window:
                read_data["end_clip_seqs"].append(left_soft_clip_seq)
        read_data["map_quals"].append(read.mapping_quality)
        read_data["strands"].append(read_strand)

    read_starts = np.array(read_data["starts"]) if read_data["starts"] else np.array([])
    read_ends = np.array(read_data["ends"]) if read_data["ends"] else np.array([])
    start_soft_clips = (
        np.array(read_data["start_soft_clips"]) if read_data["start_soft_clips"] else np.array([])
    )
    end_soft_clips = np.array(read_data["end_soft_clips"]) if read_data["end_soft_clips"] else np.array([])
    map_quals = np.array(read_data["map_quals"]) if read_data["map_quals"] else np.array([])
    strands = np.array(read_data["strands"]) if read_data["strands"] else np.array([])

    total_reads = len(valid_reads)
    strand_count = {"forward": int(np.sum(strands == strand)), "reverse": int(np.sum(strands != strand))}
    read_start_density = int(np.sum(np.abs(read_starts - pos) <= cfg.density_window)) if len(read_starts) > 0 else 0
    read_end_density = int(np.sum(np.abs(read_ends - pos) <= cfg.density_window)) if len(read_ends) > 0 else 0

    nearest_splice, softclip_bias = nearest_splice_and_softclip_bias(narrow_reads, pos, cfg, strand)
    start_entropy, end_entropy = read_start_end_entropy(read_starts, read_ends, pos, cfg)
    gradient_sharpness, max_gradient, local_cov = calculate_coverage_gradient_sharpness(
        narrow_reads, pos, cfg, strand
    )
    start_clip_sequences = read_data["start_clip_seqs"]
    end_clip_sequences = read_data["end_clip_seqs"]
    all_clips = start_clip_sequences + end_clip_sequences
    avg_clip_length = float(np.mean([len(s) for s in all_clips])) if all_clips else 0.0
    max_clip_length = float(np.max([len(s) for s in all_clips])) if all_clips else 0.0
    num_clips = float(len(all_clips))

    features: dict[str, object] = {
        "transcript_id": site.transcript_id,
        "site_type": site.site_type,
        "chrom": chrom,
        "position": pos,
        "strand": strand,
        "total_reads": total_reads,
        "read_start_density": read_start_density,
        "read_end_density": read_end_density,
        "start_soft_clip_mean": float(np.mean(start_soft_clips)) if len(start_soft_clips) > 0 else 0.0,
        "end_soft_clip_mean": float(np.mean(end_soft_clips)) if len(end_soft_clips) > 0 else 0.0,
        "start_soft_clip_max": float(np.max(start_soft_clips)) if len(start_soft_clips) > 0 else 0.0,
        "end_soft_clip_max": float(np.max(end_soft_clips)) if len(end_soft_clips) > 0 else 0.0,
        "start_soft_clip_median": float(np.median(start_soft_clips)) if len(start_soft_clips) > 0 else 0.0,
        "end_soft_clip_median": float(np.median(end_soft_clips)) if len(end_soft_clips) > 0 else 0.0,
        "start_soft_clip_count": int(len(start_soft_clips)),
        "end_soft_clip_count": int(len(end_soft_clips)),
        "mean_mapq": float(np.mean(map_quals)) if len(map_quals) > 0 else 0.0,
        "std_mapq": float(np.std(map_quals)) if len(map_quals) > 0 else 0.0,
        "strand_ratio": strand_count["forward"] / max(strand_count["reverse"], 1),
        "coverage_before": coverage_before,
        "coverage_after": coverage_after,
        "delta_coverage": delta_coverage,
        "nearest_splice_dist": nearest_splice,
        "softclip_bias": softclip_bias,
        "start_entropy": start_entropy,
        "end_entropy": end_entropy,
        "coverage_gradient_sharpness": gradient_sharpness,
        "max_coverage_gradient": max_gradient,
        "local_avg_coverage": local_cov,
        "avg_clip_length": avg_clip_length,
        "max_clip_length": max_clip_length,
        "num_clips": num_clips,
    }
    features.update(
        extract_stage1_stat_features(
            read_starts,
            start_soft_clips,
            end_soft_clips,
            coverage_before,
            coverage_after,
            total_reads,
            pos,
        )
    )
    sequence_features = softclip_kmer_features(all_clips, "", cfg.k_values)
    features.update(sequence_features)

    basic_order = [
        "transcript_id",
        "site_type",
        "chrom",
        "position",
        "strand",
        "total_reads",
        "read_start_density",
        "read_end_density",
        "start_soft_clip_mean",
        "end_soft_clip_mean",
        "start_soft_clip_max",
        "end_soft_clip_max",
        "start_soft_clip_median",
        "end_soft_clip_median",
        "start_soft_clip_count",
        "end_soft_clip_count",
        "mean_mapq",
        "std_mapq",
        "strand_ratio",
        "coverage_before",
        "coverage_after",
        "delta_coverage",
        "nearest_splice_dist",
        "softclip_bias",
        "start_entropy",
        "end_entropy",
        "coverage_gradient_sharpness",
        "max_coverage_gradient",
        "local_avg_coverage",
        "avg_clip_length",
        "max_clip_length",
        "num_clips",
    ]
    ordered = OrderedDict()
    for key in basic_order:
        if key in features:
            ordered[key] = features[key]
    for key in sorted(k for k in features.keys() if k not in ordered):
        ordered[key] = features[key]
    return dict(ordered)


def _stage1_pool_worker(task: tuple[CandidateSite, Stage1FeatureConfig]) -> dict[str, object]:
    """Picklable worker: unpack ``(site, cfg)`` and call :func:`compute_stage1_features_for_site` on pool BAM."""
    site, site_cfg = task
    if _STAGE1_POOL_BAM is None:
        raise RuntimeError("Stage I worker pool BAM handle is not initialized")
    return compute_stage1_features_for_site(_STAGE1_POOL_BAM, site, site_cfg)  # type: ignore[arg-type]


def _stage1_progress_step(n: int) -> int:
    """Print roughly every ~10% of sites, at least every 1000, at least every site when tiny."""
    if n <= 0:
        return 1
    return max(1, min(10000, n // 10))


def compute_stage1_features(
    bam_path: Path,
    candidates: Iterable[CandidateSite],
    cfg: Stage1FeatureConfig,
    *,
    parallel: bool = True,
    parallel_min_sites: int = 100,
    n_workers: int | None = None,
    progress: bool = True,
) -> list[dict[str, object]]:
    """
    Extract Stage I rows in candidate order.

    When ``parallel`` is True and there are more than ``parallel_min_sites`` candidates,
    uses a process pool with one shared BAM handle per worker.

    Set ``progress=False`` to silence ``[telos] Stage I features`` log lines.
    """
    sites = list(candidates)
    if not sites:
        return []

    bam_str = str(bam_path.resolve())
    if n_workers is None:
        n_proc = min(cpu_count(), 8)
    else:
        n_proc = max(1, int(n_workers))

    n = len(sites)
    step = _stage1_progress_step(n)
    use_pool = parallel and n > parallel_min_sites and n_proc > 1
    if not use_pool:
        if progress:
            print(f"[telos] Stage I features: {n} sites (sequential)", flush=True)
        rows: list[dict[str, object]] = []
        with pysam.AlignmentFile(bam_str, "rb") as bam:
            for i, s in enumerate(sites, 1):
                rows.append(compute_stage1_features_for_site(bam, s, cfg))
                if progress and (i == n or i % step == 0):
                    print(f"[telos] Stage I features: {i}/{n} sites", flush=True)
        if progress:
            print(f"[telos] Stage I features: finished {n} sites.", flush=True)
        return rows

    tasks = [(s, cfg) for s in sites]
    chunksize = max(1, n // (n_proc * 8))
    if progress:
        print(f"[telos] Stage I features: {n} sites ({n_proc} workers)", flush=True)
    rows_parallel: list[dict[str, object]] = []
    with Pool(processes=n_proc, initializer=_stage1_pool_init, initargs=(bam_str,)) as pool:
        for i, row in enumerate(pool.imap(_stage1_pool_worker, tasks, chunksize=chunksize), 1):
            rows_parallel.append(row)
            if progress and (i == n or i % step == 0):
                print(f"[telos] Stage I features: {i}/{n} sites", flush=True)
    if progress:
        print(f"[telos] Stage I features: finished {n} sites.", flush=True)
    return rows_parallel


def write_feature_tsv(rows: list[dict[str, object]], out_tsv: Path) -> None:
    """
    Write Stage I feature rows as TSV; header from first row keys or minimal id columns if empty.

    Uses :class:`csv.DictWriter` with tab delimiter.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        with out_tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["transcript_id", "site_type", "chrom", "position", "strand"])
        return

    fieldnames = list(rows[0].keys())
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
