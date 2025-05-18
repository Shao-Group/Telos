import numpy as np
import argparse
import pandas as pd
import os
from tqdm import tqdm
from ml_utils import normalize_chrom
from config import load_config

def load_candidates(candidate_path):
    """Load extracted candidate features."""
    return pd.read_csv(candidate_path, dtype={"chrom": str, "position": int, "strand": str})

def load_reference(reference_path):
    """Load reference annotation TSS/TES positions."""
    ref_df = pd.read_csv(reference_path, sep=" ", header=None,
                         names=["site_type", "chrom", "position", "pos_strand_cnt", "neg_strand_cnt"],
                         dtype={"chrom": str})
    ref_df["strand"] = ref_df.apply(
        lambda x: "+" if x["pos_strand_cnt"] > x["neg_strand_cnt"] else "-", axis=1
    )
    ref_df = ref_df[["site_type", "chrom", "position", "strand"]]
    # print(ref_df.head(5))
    return ref_df



def label_candidates(candidate_df, reference_df, site_type, max_distance=50):
    """Fast labeling using grouped and vectorized search."""
    
    reference_df = reference_df[
        reference_df['site_type'].str.lower() == site_type.lower()
    ].copy()
    reference_df['chrom'] = reference_df['chrom'].astype(str)
    reference_df['position'] = reference_df['position'].astype(int)

    # Group reference sites for fast lookup
    grouped_ref = reference_df.groupby(['chrom', 'strand'])

    labels = np.zeros(len(candidate_df), dtype=int)

    print(f"Labeling {site_type.upper()} candidates with {max_distance} nt tolerance ...")
    for i, row in tqdm(candidate_df.iterrows(), total=len(candidate_df), desc=f"Labeling {site_type.upper()}"):
        chrom, pos, strand = row['chrom'], row['position'], row['strand']
        
        try:
            group = grouped_ref.get_group((chrom, strand))
        except KeyError:
            continue  # No matching chrom+strand in reference

        # Use numpy for fast range check
        ref_positions = group['position'].values
        match_found = np.any(np.abs(ref_positions - pos) <= max_distance)
        labels[i] = 1 if match_found else 0

    candidate_df['label'] = labels
    return candidate_df


def main(args, cfg):
    # Load configuration

    reference_df = load_reference(cfg.ref_candidate_file)
    # Load candidates
    tss_candidate_df = load_candidates(cfg.tss_feature_file)
    tes_candidate_df = load_candidates(cfg.tes_feature_file)

    # Label candidates
    tss_labeled_df = label_candidates(tss_candidate_df, reference_df, "TSS", args.distance)
    tes_labeled_df = label_candidates(tes_candidate_df, reference_df, "TES", args.distance)

    # Save labeled candidates
    tss_labeled_df.to_csv(cfg.tss_labeled_file, index=False)
    tes_labeled_df.to_csv(cfg.tes_labeled_file, index=False)
    print(f"Labeled TSS candidates saved to: {cfg.tss_labeled_file}")
    print(f"Labeled TES candidates saved to: {cfg.tes_labeled_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Label candidate TSS/TES sites based on reference.")
    parser.add_argument('-d', '--distance', type=int, default=50, help="Maximum distance allowed for matching.")
    parser.add_argument('-c', '--config', type=str, required=True, help="Path to the configuration file.")
    args = parser.parse_args()
    cfg = load_config(args.config)
    main(args, cfg)
