import numpy as np
import argparse
import pandas as pd
import os
from tqdm import tqdm
from utils.ml_utils import normalize_chrom

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
    print(ref_df.head(5))
    return ref_df



def label_candidates(candidate_df, reference_df, site_type, max_distance=50):
    """Fast labeling using grouped and vectorized search."""
    # Ensure types are aligned
    reference_df = reference_df[
        reference_df['site_type'].str.lower() == site_type.lower()
    ].copy()
    reference_df['chrom'] = reference_df['chrom'].astype(str)
    reference_df['position'] = reference_df['position'].astype(int)

    # Group reference sites for fast lookup
    grouped_ref = reference_df.groupby(['chrom', 'strand'])

    labels = np.zeros(len(candidate_df), dtype=int)

    print(candidate_df.head(5))
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


def process_one_file(method, site_type, feature_dir, output_dir, reference_df, max_distance):
    input_path = os.path.join(feature_dir, f"{method}_{site_type}.csv")
    output_path = os.path.join(output_dir, f"{method}_{site_type}_labeled.csv")

    print(f"\nProcessing: {method.upper()} {site_type.upper()}")
    candidate_df = load_candidates(input_path)

    labeled_df = label_candidates(candidate_df, reference_df, site_type, max_distance)
    labeled_df.to_csv(output_path, index=False)
    print(f"Saved to: {output_path}")

def main(args):
    reference_df = load_reference(args.reference)
    if args.mapping:
        chrom_mapping = pd.read_csv(args.mapping, sep='\t', header=None, names=["chrom", "mapped"])
        chrom_mapping = dict(zip(chrom_mapping["chrom"], chrom_mapping["mapped"]))
        reference_df['chrom'] = reference_df['chrom'].map(chrom_mapping)
    # print reference_df after shuffling 
    # print(reference_df.sample(frac=1).head(5))
    os.makedirs(args.output_dir, exist_ok=True)

    if args.batch:
        methods = ["isoquant", "stringtie"] if args.assembly == "long" else ["stringtie", "scallop2"]
        site_types = ["tss", "tes"]
        for method in methods:
            for site_type in site_types:
                process_one_file(method, site_type, args.feature_dir, args.output_dir, reference_df, args.distance)
        print("All files processed in batch mode.")
    else:
        if not args.method or not args.site_type:
            raise ValueError("For non-batch mode, --method and --site_type must be specified.")
        process_one_file(args.method, args.site_type, args.feature_dir, args.output_dir, reference_df, args.distance)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Label candidate TSS/TES sites based on reference.")
    parser.add_argument('-f', '--feature_dir', type=str, default="features", help="Directory containing feature CSVs.")
    parser.add_argument('-o', '--output_dir', type=str, default="data_train", help="Directory to save labeled output.")
    parser.add_argument('-r', '--reference', type=str, required=True, help="Path to reference annotation file (TSV).")
    parser.add_argument('-m', '--mapping', type=str, required=False, help="Path to chrom mapping file (TSV).")
    parser.add_argument('-d', '--distance', type=int, default=50, help="Maximum distance allowed for matching.")
    parser.add_argument('-b', '--batch', action='store_true', help="Process all combinations of methods and site types.")
    parser.add_argument('-M', '--method', type=str, help="Method name (e.g., isoquant, stringtie) if not in batch mode.")
    parser.add_argument('-T', '--site_type', type=str, help="Site type (tss or tes) if not in batch mode.")
    parser.add_argument('-a', '--assembly', type=str, choices=["short", "long"], required=True, help="Assembly type (short or long).")
    args = parser.parse_args()
    main(args)
