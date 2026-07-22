"""
RefSeq-relative novelty: genome1-style transcript matching and endpoint proximity.

Classifies multi-exon transcripts in any GTF as identical to RefSeq or novel (strict intron chain +
strand + boundary match). Provides GTF/ranked-TSV filtering and endpoint-novel site flags.
"""

from __future__ import annotations

import pickle
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

import pandas as pd

from telos.backends.gtfformat import build_cov_dataframe
from telos.benchmark.stage1_metrics import safe_aupr
from telos.gtf_attributes import parse_transcript_id
from telos.labels.site_labels import (
    label_sites_by_proximity,
    normalize_chrom_name,
    reference_sites_from_gtf,
)

IntronChain = tuple[tuple[int, int], ...]
BoundaryEntry = tuple[IntronChain, str]


@dataclass
class TranscriptModel:
    transcript_id: str
    chrom: str
    strand: str
    exons: list[tuple[int, int]]


@dataclass
class RefSeqTranscriptIndex:
    """Boundary hash: (chrom_norm, first_intron_donor, last_intron_acceptor) -> chain/strand set."""

    by_boundary: dict[tuple[str, int, int], set[BoundaryEntry]] = field(default_factory=dict)
    chroms: set[str] = field(default_factory=set)

    def add(self, tx: TranscriptModel) -> None:
        chain = _intron_chain(tx.exons)
        if not chain:
            return
        chrom = normalize_chrom_name(tx.chrom)
        p1 = chain[0][0]
        p2 = chain[-1][1]
        key = (chrom, p1, p2)
        self.chroms.add(chrom)
        entry: BoundaryEntry = (chain, tx.strand)
        self.by_boundary.setdefault(key, set()).add(entry)

    def classify(self, tx: TranscriptModel) -> tuple[str, int]:
        """
        Return (match_type, is_novel_vs_refseq).

        match_type: identical | single_exon_skipped | chrm_mismatch | boundary_mismatch |
                    chain_mismatch | strand_mismatch
        """
        if len(tx.exons) < 2:
            return "single_exon_skipped", 0
        chain = _intron_chain(tx.exons)
        if not chain:
            return "single_exon_skipped", 0
        chrom = normalize_chrom_name(tx.chrom)
        if chrom not in self.chroms:
            return "chrm_mismatch", 1
        p1 = chain[0][0]
        p2 = chain[-1][1]
        key = (chrom, p1, p2)
        candidates = self.by_boundary.get(key)
        if not candidates:
            return "boundary_mismatch", 1
        saw_chain_match = False
        for ref_chain, ref_strand in candidates:
            if ref_chain != chain:
                continue
            saw_chain_match = True
            if ref_strand == tx.strand:
                return "identical", 0
        if saw_chain_match:
            return "strand_mismatch", 1
        return "chain_mismatch", 1


def _intron_chain(exons: list[tuple[int, int]]) -> IntronChain:
    if len(exons) < 2:
        return tuple()
    ordered = sorted(exons, key=lambda e: e[0])
    pairs: list[tuple[int, int]] = []
    for i in range(len(ordered) - 1):
        pairs.append((ordered[i][1], ordered[i + 1][0]))
    return tuple(pairs)


def _iter_exon_rows(gtf_path: Path) -> Iterator[tuple[str, str, str, int, int]]:
    with gtf_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            chrom, _src, feature, start, end, _score, strand, _frame, attrs = cols[:9]
            if feature != "exon":
                continue
            tid = parse_transcript_id(attrs)
            if not tid:
                continue
            yield tid, chrom, strand, int(start), int(end)


def parse_multi_exon_transcripts(gtf_path: Path) -> dict[str, TranscriptModel]:
    """Build transcript models from exon rows (multi-exon only returned)."""
    buckets: dict[str, list[tuple[str, str, int, int]]] = {}
    for tid, chrom, strand, start, end in _iter_exon_rows(gtf_path.resolve()):
        buckets.setdefault(tid, []).append((chrom, strand, start, end))
    out: dict[str, TranscriptModel] = {}
    for tid, rows in buckets.items():
        if len(rows) < 2:
            continue
        chrom = rows[0][0]
        strand_counts: dict[str, int] = {}
        exons: list[tuple[int, int]] = []
        for c, s, st, en in rows:
            chrom = c
            strand_counts[s] = strand_counts.get(s, 0) + 1
            exons.append((st, en))
        strand = max(strand_counts, key=strand_counts.get) if strand_counts else "."
        out[tid] = TranscriptModel(
            transcript_id=tid,
            chrom=chrom,
            strand=strand,
            exons=exons,
        )
    return out


def build_refseq_transcript_index(refseq_gtf: Path) -> RefSeqTranscriptIndex:
    idx = RefSeqTranscriptIndex()
    for tx in parse_multi_exon_transcripts(refseq_gtf).values():
        idx.add(tx)
    return idx


def save_refseq_index(index: RefSeqTranscriptIndex, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(index, fh, protocol=pickle.HIGHEST_PROTOCOL)


def load_refseq_index(path: Path) -> RefSeqTranscriptIndex:
    with path.open("rb") as fh:
        obj = pickle.load(fh)
    if not isinstance(obj, RefSeqTranscriptIndex):
        raise TypeError(f"Not a RefSeqTranscriptIndex: {path}")
    return obj


def get_or_build_refseq_index(refseq_gtf: Path, cache_path: Path) -> RefSeqTranscriptIndex:
    refseq_gtf = refseq_gtf.resolve()
    if cache_path.is_file() and cache_path.stat().st_mtime >= refseq_gtf.stat().st_mtime:
        print(f"[telos] refseq-novel: loading transcript index cache {cache_path}", flush=True)
        return load_refseq_index(cache_path)
    print(f"[telos] refseq-novel: building transcript index from {refseq_gtf} ...", flush=True)
    idx = build_refseq_transcript_index(refseq_gtf)
    save_refseq_index(idx, cache_path)
    print(
        f"[telos] refseq-novel: index built ({len(idx.by_boundary)} boundary keys) -> {cache_path}",
        flush=True,
    )
    return idx


def get_or_build_refseq_endpoint_sites(refseq_gtf: Path, cache_path: Path) -> pd.DataFrame:
    refseq_gtf = refseq_gtf.resolve()
    if cache_path.is_file() and cache_path.stat().st_mtime >= refseq_gtf.stat().st_mtime:
        print(f"[telos] refseq-novel: loading endpoint sites cache {cache_path}", flush=True)
        return pd.read_csv(cache_path, sep="\t")
    print(f"[telos] refseq-novel: extracting RefSeq TSS/TES from {refseq_gtf} ...", flush=True)
    sites = refseq_endpoint_sites(refseq_gtf)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    sites.to_csv(cache_path, sep="\t", index=False)
    print(f"[telos] refseq-novel: cached {len(sites)} RefSeq endpoint sites", flush=True)
    return sites


def classify_transcripts_in_gtf(gtf_path: Path, index: RefSeqTranscriptIndex) -> pd.DataFrame:
    txs = parse_multi_exon_transcripts(gtf_path)
    rows: list[dict[str, object]] = []
    for tid, tx in sorted(txs.items()):
        match_type, is_novel = index.classify(tx)
        chain = _intron_chain(tx.exons)
        rows.append(
            {
                "transcript_id": tid,
                "chrom": normalize_chrom_name(tx.chrom),
                "strand": tx.strand,
                "n_exons": len(tx.exons),
                "match_type": match_type,
                "is_novel_vs_refseq": int(is_novel),
                "first_intron_donor": chain[0][0] if chain else "",
                "last_intron_acceptor": chain[-1][1] if chain else "",
            }
        )
    if not rows:
        return pd.DataFrame(
            columns=[
                "transcript_id",
                "chrom",
                "strand",
                "n_exons",
                "match_type",
                "is_novel_vs_refseq",
                "first_intron_donor",
                "last_intron_acceptor",
            ]
        )
    return pd.DataFrame(rows)


def filter_gtf_by_transcript_ids(gtf_in: Path, keep_ids: set[str], gtf_out: Path) -> int:
    """Write GTF rows for transcript/exon/CDS features whose transcript_id is in keep_ids."""
    gtf_out.parent.mkdir(parents=True, exist_ok=True)
    n_lines = 0
    with gtf_in.open("r", encoding="utf-8", errors="replace") as inp, gtf_out.open(
        "w", encoding="utf-8"
    ) as out:
        out.write(f"## filtered_by_novel_vs_refseq n_ids={len(keep_ids)}\n")
        for raw in inp:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            feature = cols[2]
            if feature not in ("transcript", "exon", "CDS", "cds"):
                continue
            tid = parse_transcript_id(cols[8])
            if tid and tid in keep_ids:
                out.write(line + "\n")
                n_lines += 1
    return n_lines


def filter_ranked_tsv(ranked_tsv: Path, keep_ids: set[str], out_tsv: Path) -> int:
    df = pd.read_csv(ranked_tsv, sep="\t", dtype={"transcript_id": str}, low_memory=False)
    if "transcript_id" not in df.columns:
        raise ValueError(f"ranked TSV missing transcript_id: {ranked_tsv}")
    sub = df[df["transcript_id"].isin(keep_ids)].copy()
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    sub.to_csv(out_tsv, sep="\t", index=False)
    return int(len(sub))


def refseq_endpoint_sites(refseq_gtf: Path) -> pd.DataFrame:
    df = reference_sites_from_gtf(refseq_gtf.resolve())
    if df.empty:
        return df
    out = df.copy()
    out["chrom"] = out["chrom"].map(normalize_chrom_name)
    return out


def endpoint_is_novel_vs_refseq(
    sites_df: pd.DataFrame,
    refseq_sites: pd.DataFrame,
    site_type: str,
    tolerance_bp: int,
) -> pd.Series:
    """1 if site is NOT within tolerance_bp of any RefSeq endpoint of the same type."""
    matched = label_sites_by_proximity(sites_df, refseq_sites, site_type, tolerance_bp)
    return (1 - matched.astype(int)).astype(int)


def _join_sites_with_cov(site_type: str, sites_df: pd.DataFrame, assembly_gtf: Path) -> pd.DataFrame:
    cov = build_cov_dataframe(assembly_gtf)
    if cov.empty:
        return pd.DataFrame()
    st = site_type.upper()
    sub = sites_df[sites_df["site_type"].astype(str).str.upper() == st].copy()
    if sub.empty:
        return pd.DataFrame()
    sub["chrom"] = sub["chrom"].map(normalize_chrom_name)
    if st == "TSS":
        csub = cov[["tss_chrom", "tss_pos", "strand", "coverage"]].copy()
        csub["tss_chrom"] = csub["tss_chrom"].map(normalize_chrom_name)
        if len(csub) > csub[["tss_chrom", "tss_pos", "strand"]].drop_duplicates().shape[0]:
            csub = csub.groupby(["tss_chrom", "tss_pos", "strand"], as_index=False)["coverage"].mean()
        return sub.merge(
            csub,
            left_on=["chrom", "position", "strand"],
            right_on=["tss_chrom", "tss_pos", "strand"],
            how="inner",
        )
    csub = cov[["tes_chrom", "tes_pos", "strand", "coverage"]].copy()
    csub["tes_chrom"] = csub["tes_chrom"].map(normalize_chrom_name)
    if len(csub) > csub[["tes_chrom", "tes_pos", "strand"]].drop_duplicates().shape[0]:
        csub = csub.groupby(["tes_chrom", "tes_pos", "strand"], as_index=False)["coverage"].mean()
    return sub.merge(
        csub,
        left_on=["chrom", "position", "strand"],
        right_on=["tes_chrom", "tes_pos", "strand"],
        how="inner",
    )


def evaluate_novel_sites_aupr(
    *,
    sites_scored_tsv: Path,
    assembly_gtf: Path,
    test_ref_gtf: Path,
    refseq_sites: pd.DataFrame,
    tolerance_bp: int,
) -> dict[str, dict[str, float | int | None]]:
    """
    AUPR on all endpoint-novel sites (is_novel==1 vs RefSeq).

    Returns nested dict keyed by site type lower (tss/tes) then metric name.
    """
    df = pd.read_csv(
        sites_scored_tsv,
        sep="\t",
        low_memory=False,
        dtype={
            "site_type": str,
            "chrom": str,
            "position": "Int64",
            "strand": str,
            "p_site_rf": float,
            "p_site_xgb": float,
        },
    )
    need = {"site_type", "chrom", "position", "strand", "p_site_rf", "p_site_xgb"}
    if df.empty or not need.issubset(df.columns):
        return {}
    anno_sites = reference_sites_from_gtf(test_ref_gtf.resolve())
    if anno_sites.empty:
        return {}

    out: dict[str, dict[str, float | int | None]] = {}
    for st in ("TSS", "TES"):
        sub = _join_sites_with_cov(st, df, assembly_gtf)
        if sub.empty:
            continue
        labels = label_sites_by_proximity(sub, anno_sites, st, tolerance_bp).astype(int)
        is_novel = endpoint_is_novel_vs_refseq(sub, refseq_sites, st, tolerance_bp)
        nov = sub[is_novel == 1].copy()
        if nov.empty:
            continue
        nov_labels = labels.loc[nov.index]
        key = st.lower()
        n_novel = int(len(nov))
        n_pos = int(nov_labels.sum())
        out[key] = {
            "n_novel_sites": n_novel,
            "n_novel_pos": n_pos,
            "aupr_novel_rf": safe_aupr(nov_labels, nov["p_site_rf"]),
            "aupr_novel_xgb": safe_aupr(nov_labels, nov["p_site_xgb"]),
            "aupr_novel_baseline": safe_aupr(nov_labels, nov["coverage"]),
        }
    return out
