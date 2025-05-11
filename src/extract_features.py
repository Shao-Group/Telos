import argparse
import pysam
import pandas as pd
import numpy as np
from collections import Counter
from config import config
import math
import os
from tqdm import tqdm
from scipy.stats import entropy

# Example candidate sites text
# TSS GL000194.1 115066 0 1
# TSS GL000195.1 86726 0 2
# TSS GL000195.1 137958 4 0
# TSS GL000195.1 142050 1 0



def load_candidate_sites(file_path):
    with open(file_path, 'r') as f:
        candidate_sites_text = f.read()
        # Parse the candidate sites
    tss_candidate_sites = []
    tes_candidate_sites = []
    for line in candidate_sites_text.strip().split("\n"):
        parts = line.strip().split()
        site_type, chrom, pos, pos_strand_count, neg_strand_count = parts[0], parts[1], int(parts[2]), int(parts[3]), int(parts[4])
        strand = "+" if pos_strand_count > neg_strand_count else "-" 
        if site_type == "TSS":
            tss_candidate_sites.append((chrom, pos, strand))
        elif site_type == "TES":
            tes_candidate_sites.append((chrom, pos, strand))

    return tss_candidate_sites, tes_candidate_sites



def extract_features(bam, chrom, pos, strand, cfg):
    region_start = max(0, pos - cfg.window_size)
    region_end = pos + cfg.window_size
    reads = bam.fetch(chrom, region_start, region_end)

    read_starts, read_ends, soft_clips, map_quals = [], [], [], []
    strand_count = Counter()
    total_reads = 0

    for read in reads:
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
            continue
        total_reads += 1
        read_starts.append(read.reference_start)
        read_ends.append(read.reference_end)
        map_quals.append(read.mapping_quality)
        strand_count["+" if not read.is_reverse else "-"] += 1

        if read.cigartuples:
            if read.cigartuples[0][0] == 4:  # Soft clip at start
                soft_clips.append(read.cigartuples[0][1])
            if read.cigartuples[-1][0] == 4:  # Soft clip at end
                soft_clips.append(read.cigartuples[-1][1])

    window_radius = cfg.density_window
    read_start_density = sum(1 for s in read_starts if abs(s - pos) <= window_radius)
    read_end_density   = sum(1 for e in read_ends if abs(e - pos) <= window_radius)

    
    coverage_before, coverage_after, delta_coverage = calculate_coverage_change(bam, chrom, pos, cfg)
    nearest_splice = find_nearest_splice_site(bam, chrom, pos, cfg)
    softclip_bias = soft_clip_bias(bam, chrom, pos, cfg)
    start_entropy, end_entropy = read_start_end_entropy(read_starts, read_ends, pos, cfg)
    return {
        "chrom": chrom,
        "position": pos,
        "strand": strand,
        "total_reads": total_reads,
        "read_start_density": read_start_density,
        "read_end_density": read_end_density,
        "soft_clip_mean": np.mean(soft_clips) if soft_clips else 0,
        "soft_clip_max": max(soft_clips) if soft_clips else 0,
        "soft_clip_median": np.median(soft_clips) if soft_clips else 0,
        "soft_clip_count": len(soft_clips),
        "soft_clip_entropy": softclip_entropy(reads, pos, cfg.soft_clip_window),
        "mean_mapq": np.mean(map_quals) if map_quals else 0,
        "std_mapq": np.std(map_quals) if map_quals else 0,
        "strand_ratio": strand_count["+"] / max(strand_count["-"], 1),
        "coverage_before": coverage_before,
        "coverage_after": coverage_after,
        "delta_coverage": delta_coverage,
        "nearest_splice_dist": nearest_splice,
        "softclip_bias": softclip_bias,
        "start_entropy": start_entropy,
        "end_entropy": end_entropy, 
        "full_length_reads": sum(1 for r in reads if r.reference_start < region_start + 10 and r.reference_end > region_end - 10)
    }


def calculate_coverage_change(bam, chrom, pos, cfg):
    """Calculate coverage before and after the candidate site."""
    if pos < cfg.coverage_window:
        print(f"Warning: Position {pos} is too close to the start of the chromosome {chrom}.")
        coverage_before = 0
    else:
        region_before = bam.count_coverage(
            chrom, pos - cfg.coverage_window, pos, quality_threshold=cfg.min_mapq)
        coverage_before = int(np.sum(region_before))
    
    chrom_length = bam.get_reference_length(chrom)
    if pos + cfg.coverage_window > chrom_length:
        print(f"Warning: Position {pos} is too close to the end of the chromosome {chrom}.")
        coverage_after = 0
    else:
        # Ensure we don't go out of bounds
        region_after = bam.count_coverage(
            chrom, pos, pos + cfg.coverage_window, quality_threshold=cfg.min_mapq)
        coverage_after = int(np.sum(region_after))


    delta_coverage = coverage_after - coverage_before
    return coverage_before, coverage_after, delta_coverage

def find_nearest_splice_site(bam, chrom, pos, cfg):
    """Find distance to the nearest splice junction (N CIGAR operation)."""
    min_distance = cfg.splice_site_window
    for read in bam.fetch(chrom, max(0, pos - cfg.splice_site_window), pos + cfg.splice_site_window):
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq or not read.cigartuples:
            continue
        curr_pos = read.reference_start
        for op, length in read.cigartuples:
            if op == 0 or op == 7 or op == 8:  # Match/Mismatch (M/=/X)
                curr_pos += length
            elif op == 3:  # 'N' = spliced region
                splice_site_start = curr_pos
                splice_site_end = curr_pos + length
                distance = min(abs(pos - splice_site_start), abs(pos - splice_site_end))
                min_distance = min(min_distance, distance)
                curr_pos += length
            else:
                curr_pos += length
    return min_distance

def soft_clip_bias(bam, chrom, pos, cfg):
    """Ratio of reads with soft-clip starting exactly at the site."""
    soft_clip_count = 0
    total_reads = 0
    for read in bam.fetch(chrom, max(0, pos - cfg.soft_clip_window), pos + cfg.soft_clip_window):
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq or not read.cigartuples:
            continue
        total_reads += 1
        if read.cigartuples[0][0] == 4 and abs(read.reference_start - pos) <= 5:
            soft_clip_count += 1
    return soft_clip_count / total_reads if total_reads else 0


def calculate_entropy(positions):
    """Calculate entropy of read starts or ends."""
    if len(positions) == 0:
        return 0
    count = Counter(positions)
    total = sum(count.values())
    entropy = -sum((freq/total) * math.log2(freq/total) for freq in count.values())
    return entropy

def softclip_entropy(reads, pos, window=5):
    clipped_bases = []
    for r in reads:
        if not r.cigartuples: continue
        if r.cigartuples[0][0] == 4 and abs(r.reference_start - pos) <= window:
            clipped_bases.extend(r.query_sequence[:r.cigartuples[0][1]])
        if r.cigartuples[-1][0] == 4 and abs(r.reference_end - pos) <= window:
            clipped_bases.extend(r.query_sequence[-r.cigartuples[-1][1]:])
    if not clipped_bases:
        return 0
    return entropy(Counter(clipped_bases).values(), base=2)


def compute_strand(read):
    ts = read.get_tag("ts")
    if ts not in ["+", "-"]:
        # raise ValueError(f"Invalid strand tag: {ts}")
        return "."
    if read.is_reverse:
        return "-" if ts == "+" else "+"
    else:
        return "+" if ts == "+" else "-"

def read_start_end_entropy(start_positions, end_positions, pos, cfg):
    """Entropy around read start and end positions."""
    start_entropy = calculate_entropy(start_positions)
    end_entropy = calculate_entropy(end_positions)
    return start_entropy, end_entropy


def main(cfg):
    if cfg == None:
        return
    
    bam = pysam.AlignmentFile(cfg.bam_file, "rb")  # <-- adjust path if needed
    tss_candidate_sites, tes_candidate_sites = load_candidate_sites(cfg.candidate_sites_file)
    
    # Collect features
    print("Extracting TSS candidate features...")
    tss_feature_list = [extract_features(bam, *site, cfg) for site in tqdm(tss_candidate_sites, desc="TSS Feature Extraction")]
    features_df = pd.DataFrame(tss_feature_list)
    features_df.to_csv(cfg.tss_output_file, index=False)
    
    print("Extracting TES candidate features...")
    tes_feature_list = [extract_features(bam, *site, cfg) for site in tqdm(tes_candidate_sites, desc="TES Feature Extraction")]
    
    features_df = pd.DataFrame(tes_feature_list)
    features_df.to_csv(cfg.tes_output_file, index=False)
    
    # Close the BAM file
    bam.close()

    print("Feature extraction complete! Output saved as 'candidate_site_features.csv'.")


if __name__ == "__main__":
    main()