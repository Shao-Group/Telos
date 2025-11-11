import argparse
import pysam
import pandas as pd
import numpy as np
from collections import Counter, OrderedDict
from config import load_config, save_config
import math
import sys
import re
import gc
import os
from tqdm import tqdm
from scipy.stats import entropy
from sklearn.feature_selection import SelectKBest, f_classif, mutual_info_classif, RFECV
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
import gen_baseline_labels as gen_baseline_labels
# from sequence_embeddings import get_embedder, cleanup_embedder
from kmer_embeddings import get_separate_kmer_embedders, cleanup_kmer_embedders
from multiprocessing import Pool, cpu_count
from functools import partial


def load_candidate_sites(file_path):
    with open(file_path, 'r') as f:
        candidate_sites_text = f.read()
        # Parse the candidate sites
    tss_candidate_sites = []
    tes_candidate_sites = []
    for line in candidate_sites_text.strip().split("\n"):
        parts = line.strip().split()
        site_type, chrom, pos, pos_strand_count, neg_strand_count = parts[0], parts[1], int(parts[2]), int(parts[3]), int(parts[4])
        if pos_strand_count > neg_strand_count:
            strand = "+"
        elif pos_strand_count < neg_strand_count:
            strand = "-"
        else:
            strand = "."
        if site_type == "TSS":
            if strand == ".":
                tss_candidate_sites.append((chrom, pos, "+"))
                tss_candidate_sites.append((chrom, pos, "-"))
            else:
                tss_candidate_sites.append((chrom, pos, strand))
        elif site_type == "TES":
            if strand == ".":
                tes_candidate_sites.append((chrom, pos, "+"))
                tes_candidate_sites.append((chrom, pos, "-"))
            else:
                tes_candidate_sites.append((chrom, pos, strand))

    return tss_candidate_sites, tes_candidate_sites


def extract_features_parallel_worker(args):
    """Worker function for parallel processing of feature extraction."""
    site, bam_file, cfg = args
    
    # Each worker opens its own BAM file handle (required for multiprocessing)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        return extract_features(bam, *site, cfg)


def process_sites_parallel(sites, bam_file, cfg, site_type, n_processes=None):
    """Process sites in parallel using multiprocessing."""
    if n_processes is None:
        n_processes = min(cpu_count(), 8)  # Cap at 8 processes to avoid memory issues
    
    print(f"Processing {len(sites)} {site_type} sites using {n_processes} processes...")
    
    # Set current site type for k-mer embeddings
    cfg.current_site_type = site_type
    
    # Prepare arguments for worker processes
    worker_args = [(site, bam_file, cfg) for site in sites]
    
    # Use multiprocessing Pool for parallel execution
    with Pool(processes=n_processes) as pool:
        # Use tqdm for progress tracking
        if len(sites) > 100:
            from tqdm import tqdm
            results = list(tqdm(
                pool.imap(extract_features_parallel_worker, worker_args),
                total=len(sites),
                desc=f"Extracting {site_type} features"
            ))
        else:
            results = pool.map(extract_features_parallel_worker, worker_args)
    
    return results



def extract_features(bam, chrom, pos, strand, cfg):
    # Cache reads once for all feature calculations
    # Use extended window for gradient analysis (default 1000bp)
    extended_window = getattr(cfg, 'extended_window', 1000)
    region_start = max(0, pos - (cfg.window_size))
    region_end = pos + (cfg.window_size)
    cached_reads = list(bam.fetch(chrom, region_start, region_end))

    extended_cached_reads = list(bam.fetch(chrom, max(0, pos - extended_window), pos + extended_window))
    
    # Filter reads once and convert to numpy arrays for vectorized operations
    valid_reads = []
    read_data = {
        'starts': [], 'ends': [], 'start_soft_clips': [], 'end_soft_clips': [], 
        'map_quals': [], 'strands': [], 'start_clip_seqs': [], 'end_clip_seqs': []
    }
    
    for read in cached_reads:
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
            continue
            
        valid_reads.append(read)
        read_strand = compute_strand(read)
        left_soft_clip_length = 0
        right_soft_clip_length = 0
        left_soft_clip_seq = ""
        right_soft_clip_seq = ""
        
        if read.cigartuples and read.query_sequence:
            if read.cigartuples[0][0] == 4 : #and read.reference_start >= region_start:  # Soft clip at start
                left_soft_clip_length = read.cigartuples[0][1]
                left_soft_clip_seq = read.query_sequence[:left_soft_clip_length]
            if read.cigartuples[-1][0] == 4 : #and read.reference_end <= region_end:  # Soft clip at end
                right_soft_clip_length = read.cigartuples[-1][1]
                right_soft_clip_seq = read.query_sequence[-right_soft_clip_length:]

        if read_strand == "+":
            read_data['starts'].append(read.reference_start)
            read_data['ends'].append(read.reference_end)
            read_data['start_soft_clips'].append(left_soft_clip_length)
            read_data['end_soft_clips'].append(right_soft_clip_length)
            if left_soft_clip_seq != "" and abs(read.reference_start - pos) <= cfg.soft_clip_window:
                read_data['start_clip_seqs'].append(left_soft_clip_seq)
            if right_soft_clip_seq != "" and abs(read.reference_end - pos) <= cfg.soft_clip_window:
                read_data['end_clip_seqs'].append(right_soft_clip_seq)
        else:
            read_data['starts'].append(read.reference_end)
            read_data['ends'].append(read.reference_start)
            read_data['start_soft_clips'].append(right_soft_clip_length)
            read_data['end_soft_clips'].append(left_soft_clip_length)
            if right_soft_clip_seq != "" and abs(read.reference_end - pos) <= cfg.soft_clip_window:
                read_data['start_clip_seqs'].append(right_soft_clip_seq)
            if left_soft_clip_seq != "" and abs(read.reference_start - pos) <= cfg.soft_clip_window:
                read_data['end_clip_seqs'].append(left_soft_clip_seq)

        read_data['map_quals'].append(read.mapping_quality)
        read_data['strands'].append(read_strand)

    # Convert to numpy arrays for vectorized operations
    read_starts = np.array(read_data['starts']) if read_data['starts'] else np.array([])
    read_ends = np.array(read_data['ends']) if read_data['ends'] else np.array([])
    start_soft_clips = np.array(read_data['start_soft_clips']) if read_data['start_soft_clips'] else np.array([])
    end_soft_clips = np.array(read_data['end_soft_clips']) if read_data['end_soft_clips'] else np.array([])
    map_quals = np.array(read_data['map_quals']) if read_data['map_quals'] else np.array([])
    strands = np.array(read_data['strands']) if read_data['strands'] else np.array([])
    
    total_reads = len(valid_reads)
    strand_count = {"forward": np.sum(strands == strand), "reverse": np.sum(strands != strand)}

    # Vectorized density calculations
    window_radius = cfg.density_window
    read_start_density = np.sum(np.abs(read_starts - pos) <= window_radius) if len(read_starts) > 0 else 0
    read_end_density = np.sum(np.abs(read_ends - pos) <= window_radius) if len(read_ends) > 0 else 0

    # Use cached reads for all feature calculations
    coverage_before, coverage_after, delta_coverage = calculate_coverage_change_cached(extended_cached_reads, chrom, pos, cfg, strand)
    nearest_splice = find_nearest_splice_site_cached(cached_reads, chrom, pos, cfg)
    soft_clip_bias_value = soft_clip_bias_cached(cached_reads, chrom, pos, cfg, strand)
    start_entropy, end_entropy = read_start_end_entropy(read_starts, read_ends, pos, cfg)
    
    # Calculate coverage gradient sharpness
    gradient_sharpness, max_gradient, local_coverage = calculate_coverage_gradient_sharpness(cached_reads, chrom, pos, cfg, strand)
    
    # Use soft-clipped sequences collected during main loop
    start_clip_sequences = read_data['start_clip_seqs']
    end_clip_sequences = read_data['end_clip_seqs']
    all_clips = start_clip_sequences + end_clip_sequences
    
    # Base features
    features = {
        "chrom": chrom,
        "position": pos,
        "strand": strand,
        "total_reads": total_reads,
        "read_start_density": read_start_density,
        "read_end_density": read_end_density,
        "start_soft_clip_mean": np.mean(start_soft_clips) if len(start_soft_clips) > 0 else 0,
        "end_soft_clip_mean": np.mean(end_soft_clips) if len(end_soft_clips) > 0 else 0,
        "start_soft_clip_max": max(start_soft_clips) if len(start_soft_clips) > 0 else 0,
        "end_soft_clip_max": max(end_soft_clips) if len(end_soft_clips) > 0 else 0,
        "start_soft_clip_median": np.median(start_soft_clips) if len(start_soft_clips) > 0 else 0,
        "end_soft_clip_median": np.median(end_soft_clips) if len(end_soft_clips) > 0 else 0,
        "start_soft_clip_count": len(start_soft_clips),
        "end_soft_clip_count": len(end_soft_clips),
        "mean_mapq": np.mean(map_quals) if len(map_quals) > 0 else 0,
        "std_mapq": np.std(map_quals) if len(map_quals) > 0 else 0,
        "strand_ratio": strand_count["forward"] / max(strand_count["reverse"], 1),
        "coverage_before": coverage_before,
        "coverage_after": coverage_after,
        "delta_coverage": delta_coverage,
        "nearest_splice_dist": nearest_splice,
        "softclip_bias": soft_clip_bias_value,
        "start_entropy": start_entropy,
        "end_entropy": end_entropy,
        "coverage_gradient_sharpness": gradient_sharpness,
        "max_coverage_gradient": max_gradient,
        "local_avg_coverage": local_coverage
        # "full_length_reads": len([r for r in cached_reads if r.reference_start < region_start + 10 and r.reference_end > region_end - 10])
    }
    
    # Add dRNA-specific features
    features.update(extract_drna_specific_features(read_starts, read_ends, start_soft_clips, end_soft_clips, 
                                                  coverage_before, coverage_after, total_reads, pos, cfg))
    
    # Add sequence-based features (embeddings or handcrafted)
    if hasattr(cfg, 'use_embeddings') and cfg.use_embeddings:
        if hasattr(cfg, 'embedding_mode') and cfg.embedding_mode == 'hybrid':
            # Use hybrid approach: embeddings + key handcrafted features
            sequence_features = extract_hybrid_features(start_clip_sequences, end_clip_sequences, cfg)
        else:
            # Use pure embeddings approach
            raise NotImplementedError("Pure embeddings approach is not implemented yet")
            # sequence_features = extract_sequence_embeddings(start_clip_sequences, end_clip_sequences, cfg)
        features.update(sequence_features)
    else:
        # Use comprehensive handcrafted features with combined analysis only (better performance)
        # Determine site type from config
        current_site_type = getattr(cfg, 'current_site_type', 'TSS')
        combined_kmer_features = softclip_kmer_features(
            start_clip_sequences + end_clip_sequences, '', site_type=current_site_type
        )
        features.update(combined_kmer_features)
    
    # Ensure absolute order consistency by sorting feature names
    # Keep basic features first, then sort pattern-based features
    basic_features = ['chrom', 'position', 'strand', 'total_reads', 'read_start_density', 
                     'read_end_density', 'start_soft_clip_mean', 'end_soft_clip_mean',
                     'start_soft_clip_max', 'end_soft_clip_max', 'start_soft_clip_median',
                     'end_soft_clip_median', 'start_soft_clip_count', 'end_soft_clip_count',
                     'mean_mapq', 'std_mapq', 'strand_ratio', 'coverage_before',
                     'coverage_after', 'delta_coverage', 'nearest_splice_dist',
                     'softclip_bias', 'start_entropy', 'end_entropy', 'coverage_gradient_sharpness',
                     'max_coverage_gradient', 'local_avg_coverage']
    
    ordered_features = OrderedDict()
    # Add basic features in predefined order
    for key in basic_features:
        if key in features:
            ordered_features[key] = features[key]
    
    # Add remaining features in alphabetical order for consistency
    remaining_keys = sorted([k for k in features.keys() if k not in basic_features])
    for key in remaining_keys:
        ordered_features[key] = features[key]
    
    return dict(ordered_features)


def calculate_coverage_change_cached(cached_reads, chrom, pos, cfg, strand):
    """Calculate coverage before and after the candidate site using cached reads."""
    coverage_before = 0
    coverage_after = 0
    
    for read in cached_reads:
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
            continue
            
        # Count coverage in before window
        if pos - cfg.coverage_window <= read.reference_start < pos or pos - cfg.coverage_window < read.reference_end <= pos:
            overlap_start = max(read.reference_start, pos - cfg.coverage_window)
            overlap_end = min(read.reference_end, pos)
            if overlap_end > overlap_start:
                coverage_before += overlap_end - overlap_start
                
        # Count coverage in after window  
        if pos <= read.reference_start < pos + cfg.coverage_window or pos < read.reference_end <= pos + cfg.coverage_window:
            overlap_start = max(read.reference_start, pos)
            overlap_end = min(read.reference_end, pos + cfg.coverage_window)
            if overlap_end > overlap_start:
                coverage_after += overlap_end - overlap_start

    delta_coverage = coverage_after - coverage_before
    if strand == "+":
        return coverage_before, coverage_after, delta_coverage
    else:
        return coverage_after, coverage_before, -delta_coverage

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

def find_nearest_splice_site_cached(cached_reads, chrom, pos, cfg):
    """Find distance to the nearest splice junction using cached reads."""
    min_distance = cfg.splice_site_window
    for read in cached_reads:
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

def soft_clip_bias(bam, chrom, pos, cfg, strand):
    """Ratio of reads with soft-clip starting exactly at the site."""
    soft_clip_count = 0
    total_reads = 0
    for read in bam.fetch(chrom, max(0, pos - cfg.soft_clip_window), pos + cfg.soft_clip_window):
        read_strand = compute_strand(read)
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq or not read.cigartuples or read_strand != strand:
            continue
        total_reads += 1
        if read.cigartuples[0][0] == 4 and abs(read.reference_start - pos) <= 5:
            soft_clip_count += 1
        elif read.cigartuples[-1][0] == 4 and abs(read.reference_end - pos) <= 5:
            soft_clip_count += 1
    return soft_clip_count / total_reads if total_reads else 0

def soft_clip_bias_cached(cached_reads, chrom, pos, cfg, strand):
    """Ratio of reads with soft-clip starting exactly at the site using cached reads."""
    soft_clip_count = 0
    total_reads = 0
    for read in cached_reads:
        read_strand = compute_strand(read)
        if read.is_unmapped or read.mapping_quality < cfg.min_mapq or not read.cigartuples or read_strand != strand:
            continue
        total_reads += 1
        if read.cigartuples[0][0] == 4 and abs(read.reference_start - pos) <= 5:
            soft_clip_count += 1
        elif read.cigartuples[-1][0] == 4 and abs(read.reference_end - pos) <= 5:
            soft_clip_count += 1
    return soft_clip_count / total_reads if total_reads else 0


def calculate_coverage_gradient_sharpness(cached_reads, chrom, pos, cfg, strand):
    """Calculate coverage gradient sharpness to detect sharp coverage changes."""
    # Parameters
    window_size = 10  # Small window for gradient calculation
    analysis_range = getattr(cfg, 'gradient_analysis_range', 100)  # Range around site to analyze
    
    # Calculate coverage in small windows around the site
    positions = range(pos - analysis_range, pos + analysis_range + 1, window_size)
    coverages = []
    
    for window_start in positions:
        window_end = window_start + window_size
        coverage = 0
        
        for read in cached_reads:
            if read.is_unmapped or read.mapping_quality < cfg.min_mapq:
                continue
                
            # Calculate overlap between read and window
            overlap_start = max(read.reference_start, window_start)
            overlap_end = min(read.reference_end, window_end)
            
            if overlap_end > overlap_start:
                coverage += overlap_end - overlap_start
        
        # Normalize by window size
        coverages.append(coverage / window_size)
    
    if len(coverages) < 3:
        return 0.0, 0.0, 0.0
    
    # Convert to numpy array for vectorized operations
    coverages = np.array(coverages)
    
    # Calculate gradients (differences between adjacent windows)
    gradients = np.diff(coverages)
    
    # Find maximum absolute gradient (steepest change)
    max_gradient = np.max(np.abs(gradients)) if len(gradients) > 0 else 0.0
    
    # Calculate local average coverage for normalization
    local_avg_coverage = np.mean(coverages) if len(coverages) > 0 else 1.0
    
    # Calculate sharpness: normalize gradient by local coverage
    gradient_sharpness = max_gradient / max(local_avg_coverage, 1.0)
    
    # For strand-specific interpretation
    if strand == "-":
        # For negative strand, we expect opposite gradient pattern
        # Find the gradient with correct sign for the strand
        signed_gradients = gradients if strand == "+" else -gradients
        max_signed_gradient = np.max(signed_gradients) if len(signed_gradients) > 0 else 0.0
        
        # Use signed gradient for strand-specific analysis
        if max_signed_gradient > 0:
            gradient_sharpness = max_signed_gradient / max(local_avg_coverage, 1.0)
    
    return gradient_sharpness, max_gradient, local_avg_coverage


def calculate_entropy(positions):
    """Calculate entropy of read starts or ends."""
    if len(positions) == 0:
        return 0
    count = Counter(positions)
    total = sum(count.values())
    entropy = -sum((freq/total) * math.log2(freq/total) for freq in count.values())
    return entropy



def softclip_kmer_features(clipped_sequences, prefix="", klist=[3], site_type="TSS"):
    """Optimized comprehensive soft-clip feature extraction with single-pass processing."""
    if not clipped_sequences:
        clipped_sequences = [""]
    combined = {}
    
    # Pre-process all sequences once
    processed_sequences = [seq.upper() for seq in clipped_sequences if seq]
    all_bases = ''.join(processed_sequences)
    total_bases = len(all_bases)
    
    # 1. Nucleotide composition (vectorized counting)
    if total_bases > 0:
        base_counts = np.array([all_bases.count(base) for base in 'ATGC'])
        a_count, t_count, g_count, c_count = base_counts
        
        combined[f'{prefix}gc_content'] = (g_count + c_count) / total_bases
        combined[f'{prefix}purine_ratio'] = (a_count + g_count) / total_bases
    else:
        combined[f'{prefix}gc_content'] = 0
        combined[f'{prefix}purine_ratio'] = 0
    
    # 2. Homopolymer features (optimized regex compilation)
    base_patterns = {base: re.compile(f'{base}+') for base in 'ACGT'}
    for base in 'ACGT':
        max_run = 0
        pattern = base_patterns[base]
        for seq in processed_sequences:
            matches = pattern.findall(seq)
            if matches:
                seq_max = max(len(match) for match in matches)
                max_run = max(max_run, seq_max)
        combined[f'{prefix}max_poly{base}'] = max_run
    
    # 3. Optimized K-mer features (single pass for all k-values)
    combined.update(_extract_optimized_kmer_features(processed_sequences, prefix, klist))
    
    return combined



def _extract_optimized_kmer_features(processed_sequences, prefix, klist):
    """Optimized k-mer feature extraction for multiple k-values in single pass."""
    features = {}
    
    # Fixed k-mers for biological relevance (only 3-mers)
    fixed_kmers_3 = set(['ATG', 'TAA', 'TGA', 'TAG', 'GCC', 'CGC', 'GCG', 
                        'AAA', 'TTT', 'AAT', 'GTG', 'CAG', 'GAG', 
                        'CCC', 'GGG', 'CCG'])
    
    # Process all k-values simultaneously
    kmer_data = {k: {'counts': Counter(), 'fixed_counts': Counter(), 
                     'composition_counts': {'at_rich_kmers': 0, 'balanced_kmers': 0, 
                                          'gc_rich_kmers': 0, 'purine_rich_kmers': 0, 
                                          'repeat_kmers': 0}, 
                     'total_kmers': 0, 'gc_kmers': 0} for k in klist}
    
    # Single pass through all sequences for all k-values
    for seq in processed_sequences:
        seq_len = len(seq)
        for k in klist:
            if seq_len < k:
                continue
                
            data = kmer_data[k]
            for i in range(seq_len - k + 1):
                kmer = seq[i:i+k]
                if all(base in 'ATGC' for base in kmer):
                    data['counts'][kmer] += 1
                    data['total_kmers'] += 1
                    
                    # Fixed k-mer counting (only for 3-mers)
                    if k == 3 and kmer in fixed_kmers_3:
                        data['fixed_counts'][kmer] += 1
                    
                    # Pre-calculate counts for efficiency
                    gc_count = kmer.count('G') + kmer.count('C')
                    a_count = kmer.count('A')
                    g_count = kmer.count('G')
                    
                    # GC content analysis
                    gc_content = gc_count / k
                    if gc_content >= 0.5:
                        data['gc_kmers'] += 1
                    
                    # Composition analysis
                    purine_content = (a_count + g_count) / k
                    if len(set(kmer)) == 1:
                        data['composition_counts']['repeat_kmers'] += 1
                    elif gc_content > 0.6:
                        data['composition_counts']['gc_rich_kmers'] += 1
                    elif gc_content < 0.4:
                        data['composition_counts']['at_rich_kmers'] += 1
                    else:
                        data['composition_counts']['balanced_kmers'] += 1
                    
                    if purine_content > 0.6:
                        data['composition_counts']['purine_rich_kmers'] += 1
    
    # Compile features for all k-values
    for k in klist:
        data = kmer_data[k]
        
        # Fixed k-mer features (only for 3-mers)
        if k == 3:
            for kmer in fixed_kmers_3:
                features[f'{prefix}kmer_{kmer}'] = data['fixed_counts'].get(kmer, 0)
        
        # Statistical features
        total_kmers = data['total_kmers']
        if total_kmers > 0:
            kmer_counts = data['counts']
            unique_kmers = len(kmer_counts)
            diversity = unique_kmers / total_kmers
            most_frequent_count = max(kmer_counts.values()) if kmer_counts else 0
            most_frequent_ratio = most_frequent_count / total_kmers
            
            # Shannon entropy (vectorized)
            counts_array = np.array(list(kmer_counts.values()))
            probs = counts_array / total_kmers
            kmer_entropy = -np.sum(probs * np.log2(probs))
            
            features.update({
                f'{prefix}k{k}_gc_kmers_ratio': data['gc_kmers'] / total_kmers,
                f'{prefix}k{k}_kmer_diversity': diversity,
                f'{prefix}k{k}_kmer_entropy': kmer_entropy,
                f'{prefix}k{k}_most_frequent_kmer_ratio': most_frequent_ratio,
                f'{prefix}k{k}_unique_kmers': unique_kmers
            })
            
            # Composition features
            for comp_type, count in data['composition_counts'].items():
                features[f'{prefix}k{k}_{comp_type}'] = count
        else:
            # Zero features when no k-mers found
            zero_features = [f'{prefix}k{k}_gc_kmers_ratio', f'{prefix}k{k}_kmer_diversity', 
                            f'{prefix}k{k}_kmer_entropy', f'{prefix}k{k}_most_frequent_kmer_ratio',
                            f'{prefix}k{k}_unique_kmers']
            for comp_type in data['composition_counts'].keys():
                zero_features.append(f'{prefix}k{k}_{comp_type}')
            features.update({feat: 0 for feat in zero_features})
    
    return features

def _extract_kmer_features_optimized(clipped_sequences, prefix, k):
    """Legacy function - use _extract_optimized_kmer_features instead."""
    return _extract_optimized_kmer_features(clipped_sequences, prefix, [k])

def _extract_enhanced_kmer_features(clipped_sequences, prefix, k):
    """Enhanced k-mer feature extraction with positional bias."""
    # Get standard k-mer features first
    standard_features = _extract_kmer_features_optimized(clipped_sequences, prefix, k)
    return standard_features
    


def extract_kmer_embeddings(start_clip_sequences, end_clip_sequences, cfg, site_type):
    """Extract k-mer embeddings from soft-clipped sequences."""
    try:
        # Get the appropriate k-mer embedders
        start_tss, start_tes, end_tss, end_tes = get_separate_kmer_embedders()
        
        # Determine which embedders to use based on site type
        if site_type.upper() == 'TSS':
            start_embedder = start_tss
            end_embedder = end_tss
        elif site_type.upper() == 'TES':
            start_embedder = start_tes
            end_embedder = end_tes
        else:
            # Default to TSS embedders
            start_embedder = start_tss
            end_embedder = end_tss
            print(f"Warning: Unknown site_type '{site_type}', using TSS embedders")
        
        # Generate embeddings for start and end clips
        start_embedding = start_embedder.embed_site(start_clip_sequences)
        end_embedding = end_embedder.embed_site(end_clip_sequences)
        
        # Create feature dictionary
        features = {}
        
        # Add start k-mer embeddings
        start_feature_names = start_embedder.get_feature_names()
        for i, (name, val) in enumerate(zip(start_feature_names, start_embedding)):
            features[f'start_kmer_{name}'] = float(val)
        
        # Add end k-mer embeddings  
        end_feature_names = end_embedder.get_feature_names()
        for i, (name, val) in enumerate(zip(end_feature_names, end_embedding)):
            features[f'end_kmer_{name}'] = float(val)
        
        # Add combined statistics
        combined_embedding = (start_embedding + end_embedding) / 2
        for i, val in enumerate(combined_embedding):
            features[f'combined_kmer_{i}'] = float(val)
        
        return features
        
    except Exception as e:
        print(f"Warning: Could not extract k-mer embeddings: {e}")
        # Return empty features on error
        return {}



def extract_hybrid_features(start_clip_sequences, end_clip_sequences, cfg):
    """Extract both embeddings and key handcrafted features."""
    # Get embeddings (could be CNN or k-mer based on cfg.embedding_type)
    # embedding_features = extract_kmer_embeddings(start_clip_sequences, end_clip_sequences, cfg, cfg.current_site_type)
    # Use only combined features for better performance
    current_site_type = getattr(cfg, 'current_site_type', 'TSS')
    combined_kmer_features = softclip_kmer_features(
        start_clip_sequences + end_clip_sequences, '', site_type=current_site_type
    )

    embedding_features = combined_kmer_features
    # Keep only the most important handcrafted features
    all_clips = start_clip_sequences + end_clip_sequences
    
    # Basic sequence statistics
    total_clip_length = sum(len(seq) for seq in all_clips)
    avg_clip_length = np.mean([len(seq) for seq in all_clips]) if all_clips else 0
    max_clip_length = max(len(seq) for seq in all_clips) if all_clips else 0
    num_clips = len(all_clips)
    
    # Basic composition
    if all_clips:
        all_bases = ''.join(all_clips).upper()
        total_bases = len(all_bases)
        gc_content = (all_bases.count('G') + all_bases.count('C')) / total_bases if total_bases > 0 else 0
        at_content = (all_bases.count('A') + all_bases.count('T')) / total_bases if total_bases > 0 else 0
    else:
        gc_content = 0
        at_content = 0
    
    # Key handcrafted features
    basic_features = {
        'total_clip_length': total_clip_length,
        'avg_clip_length': avg_clip_length,
        'max_clip_length': max_clip_length,
        'num_clips': num_clips,
        'gc_content_basic': gc_content,
        'at_content_basic': at_content,
    }
    
    # Combine embeddings with basic features
    features = {}
    features.update(embedding_features)
    features.update(basic_features)
    
    return features


def compute_strand(read):
    try:
        ts = read.get_tag("ts")
        if ts not in ["+", "-"]:
            return "."
        if read.is_reverse:
            return "-" if ts == "+" else "+"
        else:
            return "+" if ts == "+" else "-"
    except KeyError:
        # No ts tag present, infer from read orientation
        return "-" if read.is_reverse else "+"

def read_start_end_entropy(start_positions, end_positions, pos, cfg):
    """Entropy around read start and end positions."""
    start_entropy = calculate_entropy(start_positions)
    end_entropy = calculate_entropy(end_positions)
    return start_entropy, end_entropy


def extract_drna_specific_features(read_starts, read_ends, start_soft_clips, end_soft_clips, 
                                  coverage_before, coverage_after, total_reads, pos, cfg):
    """Extract dRNA-specific features for better TSS detection using vectorized operations."""
    features = {}
    
    # 1. Coverage normalization features (dRNA has higher coverage but different patterns)
    features['coverage_ratio'] = coverage_after / max(coverage_before, 1)
    features['coverage_log_ratio'] = math.log2(max(coverage_after, 1) / max(coverage_before, 1)) if coverage_before > 0 else 0
    features['normalized_delta_coverage'] = (coverage_after - coverage_before) / max(coverage_after + coverage_before, 1)
    
    # 2. Read start precision (dRNA has more variable start positions due to degradation)
    if len(read_starts) > 0:
        read_start_variance = np.var(read_starts)
        read_start_mean = np.mean(read_starts)
        features['read_start_variance'] = read_start_variance
        features['read_start_coefficient_variation'] = read_start_variance / max(read_start_mean, 1)
        
        # Distance-weighted read start density (vectorized)
        distances = np.abs(read_starts - pos)
        weighted_density = np.sum(1 / np.maximum(distances, 1))
        features['weighted_read_start_density'] = weighted_density
        
        # Read start clustering (vectorized)
        close_reads = np.sum(distances <= 5)
        features['read_start_clustering'] = close_reads / max(len(read_starts), 1)
    else:
        features['read_start_variance'] = 0
        features['read_start_coefficient_variation'] = 0
        features['weighted_read_start_density'] = 0
        features['read_start_clustering'] = 0
    
    # 3. Soft-clip sparsity and quality features (vectorized operations)
    total_soft_clips = len(start_soft_clips) + len(end_soft_clips)
    features['softclip_sparsity'] = total_soft_clips / max(total_reads, 1)
    
    if total_soft_clips > 0:
        # Combine soft-clip arrays for vectorized operations
        all_soft_clips = np.concatenate([start_soft_clips, end_soft_clips])
        features['softclip_length_variance'] = np.var(all_soft_clips)
        features['softclip_length_skewness'] = calculate_skewness(all_soft_clips)
        
        # Vectorized ratio calculation - all_soft_clips already contains the lengths
        long_clips = np.sum(all_soft_clips > 20)
        short_clips = np.sum(all_soft_clips <= 10)
        features['long_to_short_clip_ratio'] = long_clips / max(short_clips, 1)
    else:
        features['softclip_length_variance'] = 0
        features['softclip_length_skewness'] = 0
        features['long_to_short_clip_ratio'] = 0
    
    # 4. 5' degradation indicators (vectorized operations)
    if len(read_starts) > 0:
        # Vectorized upstream/downstream calculation
        upstream_starts = np.sum(read_starts < pos)
        downstream_starts = np.sum(read_starts > pos)
        features['upstream_downstream_ratio'] = upstream_starts / max(downstream_starts, 1)
        
        # 5' degradation gradient (reads should accumulate at true TSS)
        degradation_score = calculate_degradation_score(read_starts, pos)
        features['five_prime_degradation_score'] = degradation_score
    else:
        features['upstream_downstream_ratio'] = 0
        features['five_prime_degradation_score'] = 0
    
    # 5. Coverage-normalized features (account for dRNA's higher coverage)
    features['normalized_read_density'] = len(read_starts) / max(coverage_after, 1)
    features['coverage_efficiency'] = total_reads / max(coverage_before + coverage_after, 1)
    
    return features


def calculate_skewness(values):
    """Calculate skewness of a distribution."""
    if len(values) < 3:
        return 0
    mean_val = np.mean(values)
    std_val = np.std(values)
    if std_val == 0:
        return 0
    
    skew = np.mean([(x - mean_val)**3 for x in values]) / (std_val**3)
    return skew


def calculate_degradation_score(read_starts, pos):
    """Calculate 5' degradation score - higher score indicates more degradation."""
    if len(read_starts) == 0:
        return 0
    
    # Count reads in windows around the position
    window_size = 10
    windows = []
    for i in range(-50, 51, window_size):
        window_start = pos + i
        window_end = pos + i + window_size
        count = len([s for s in read_starts if window_start <= s < window_end])
        windows.append(count)
    
    if not windows or sum(windows) == 0:
        return 0
    
    # Calculate the slope of read density (negative slope indicates degradation)
    x = list(range(len(windows)))
    y = windows
    
    # Simple linear regression slope
    n = len(x)
    slope = (n * sum(x[i] * y[i] for i in range(n)) - sum(x) * sum(y)) / max(n * sum(x[i]**2 for i in range(n)) - sum(x)**2, 1)
    
    return -slope  # Return negative slope as degradation score





def main(cfg, config_path):
    if cfg == None:
        raise ValueError("Configuration has not been created. Please call create_config() first.")
        return

    cov_file, candidate_file, ref_candidate_file = gen_baseline_labels.main(cfg)
    cfg.set_cov_file(cov_file)
    cfg.set_candidate_file(candidate_file)
    cfg.set_ref_candidate_file(ref_candidate_file)

    bam = pysam.AlignmentFile(cfg.bam_file, "rb")  # <-- adjust path if needed
    tss_candidate_sites, tes_candidate_sites = load_candidate_sites(cfg.candidate_file)
    
    # # Collect features
    # if os.path.exists(cfg.tss_feature_file) and os.path.exists(cfg.tes_feature_file):
    #     print(f"Feature files already exist: {cfg.tss_feature_file}, {cfg.tes_feature_file}")
    #     # Close the BAM file
    #     bam.close()
    #     save_config(config_path)

    #     return

    print("Extracting TSS candidate features...")
    
    # Check if parallel processing should be used
    use_parallel = getattr(cfg, 'use_parallel_processing', True)
    n_processes = getattr(cfg, 'n_processes', None)
    
    if use_parallel and len(tss_candidate_sites) > 50:  # Only use parallel for larger datasets
        # Parallel processing
        all_features = process_sites_parallel(tss_candidate_sites, cfg.bam_file, cfg, 'TSS', n_processes)
    else:
        # Sequential processing (fallback)
        cfg.current_site_type = 'TSS'
        total_sites = len(tss_candidate_sites)
        print(f"Processing {total_sites} TSS candidate sites sequentially...")
        
        all_features = []
        for idx, site in enumerate(tss_candidate_sites):
            all_features.append(extract_features(bam, *site, cfg))
            if idx % 500 == 0:
                print(f"Extracted {idx} TSS candidate features")
                sys.stdout.flush()
    
    # Save all features to CSV in one operation
    tss_df = pd.DataFrame(all_features)
    tss_df.to_csv(cfg.tss_feature_file, index=False)
    
    # Clear memory
    del all_features
    del tss_df
    gc.collect()

    print("Extracting TES candidate features...")
    
    if use_parallel and len(tes_candidate_sites) > 50:  # Only use parallel for larger datasets
        # Parallel processing
        all_features = process_sites_parallel(tes_candidate_sites, cfg.bam_file, cfg, 'TES', n_processes)
    else:
        # Sequential processing (fallback)
        cfg.current_site_type = 'TES'
        total_sites = len(tes_candidate_sites)
        print(f"Processing {total_sites} TES candidate sites sequentially...")
        
        all_features = []
        for idx, site in enumerate(tes_candidate_sites):
            all_features.append(extract_features(bam, *site, cfg))
            if idx % 500 == 0:
                print(f"Extracted {idx} TES candidate features")
                sys.stdout.flush()
    
    # Save all features to CSV in one operation
    tes_df = pd.DataFrame(all_features)
    tes_df.to_csv(cfg.tes_feature_file, index=False)
    
    # Clear memory
    del all_features
    del tes_df
    gc.collect()
    
    # Clean up embedder to free GPU memory
    if hasattr(cfg, 'use_embeddings') and cfg.use_embeddings:
        embedding_type = getattr(cfg, 'embedding_type', 'kmer')
        if embedding_type == 'kmer':
            print("Cleaning up k-mer embedders...")
            cleanup_kmer_embedders()
    
    print("Feature extraction complete!")
    
    
    
    save_config(config_path)
    
    print("Feature extraction and selection complete!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract features from BAM file.")
    parser.add_argument(
        "--config",
        type=str,
        required=True,
        help="Path to the configuration file."
    )
    parser.add_argument(
        "--no-parallel",
        action="store_true",
        help="Disable parallel processing (use sequential processing)."
    )
    parser.add_argument(
        "--n-processes",
        type=int,
        default=None,
        help="Number of processes for parallel processing (default: auto-detect, max 8)."
    )
    args = parser.parse_args()
    cfg = load_config(args.config)
    
    # Set parallel processing options
    cfg.use_parallel_processing = not args.no_parallel
    if args.n_processes:
        cfg.n_processes = args.n_processes
    
    main(cfg, args.config)