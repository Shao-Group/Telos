import pysam
import pandas as pd
import numpy as np
from collections import Counter
import config
# Your candidate sites (based on the format you shared)
# candidate_sites_text = """
# TSS GL000194.1 115066 0 1
# TSS GL000195.1 86726 0 2
# TSS GL000195.1 137958 4 0
# TSS GL000195.1 142050 1 0
# """
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



def extract_features(bam, chrom, pos, strand, site_type, cfg):
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

    return {
        "chrom": chrom,
        "position": pos,
        "strand": strand,
        "site_type": site_type,
        "total_reads": total_reads,
        "read_start_density": read_starts.count(pos),
        "read_end_density": read_ends.count(pos),
        "soft_clip_mean": np.mean(soft_clips) if soft_clips else 0,
        "soft_clip_max": max(soft_clips) if soft_clips else 0,
        "mean_mapq": np.mean(map_quals) if map_quals else 0,
        "std_mapq": np.std(map_quals) if map_quals else 0,
        "strand_ratio": strand_count["+"] / max(strand_count["-"], 1),
    }



def main():
    cfg = config.config()
    # Open your BAM file
    bam = pysam.AlignmentFile(cfg.bam_file, "rb")  # <-- adjust path if needed
    # Load candidate sites
    tss_candidate_sites, tes_candidate_sites = load_candidate_sites(cfg.candidate_sites_file)
    # Collect features
    tss_feature_list = [extract_features(bam, *site, cfg) for site in tss_candidate_sites]
    features_df = pd.DataFrame(tss_feature_list)
    features_df.to_csv(cfg.tss_output_file, index=False)
    # Collect features for TES
    tes_feature_list = [extract_features(bam, *site, cfg) for site in tes_candidate_sites]
    features_df = pd.DataFrame(tes_feature_list)
    features_df.to_csv(cfg.tes_output_file, index=False)
    # Close the BAM file
    bam.close()
    
    print("Feature extraction complete! Output saved as 'candidate_site_features.csv'.")


if __name__ == "__main__":
    main()