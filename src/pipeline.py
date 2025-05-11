from config import config
import os
import argparse
import extract_features
import subprocess
import sys
import train_all
import utils.generate_roc_data as generate_roc_data
from utils.generate_pr_curves import plot_pr_curves

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Extract features from BAM file.")
    parser.add_argument("-m","--method", type=str, required=True, choices=["short", "long"], help="Method to use for candidate site extraction.")
    parser.add_argument("-b","--bam_file", type=str, required=True, help="Path to the BAM file.")
    parser.add_argument("-c","--candidate_sites_folder", type=str, help="Path to the candidate sites folder.")
    parser.add_argument("-d","--data_name", type=str, required=True, help="Name of the data set.")
    args = parser.parse_args()

    if args.candidate_sites_folder is None:
        args.candidate_sites_folder = os.path.join("data", args.data_name)

    # Initialize configuration
    short_read_assemblers = ["stringtie", "scallop2"]
    long_read_assemblers = ["stringtie","isoquant"]
    assemblers = short_read_assemblers if args.method == "short" else long_read_assemblers
    
    reference_file = "data/refSeq.tsstes"
    chrom_mapping = "data/GRCh38_RefSeq2UCSC.txt"
    feature_dir = f"features/{args.data_name}"
    train_dir = f"data_train/{args.data_name}"
    os.makedirs(feature_dir, exist_ok=True)
    os.makedirs(train_dir, exist_ok=True)
    

    print (f"🔍Feature Extraction: {assemblers} ---> {args.bam_file} ---> {args.candidate_sites_folder} ---> {args.data_name}")
    for assembler in assemblers:
        candidate_sites_file = os.path.join(args.candidate_sites_folder, f"{args.data_name}_{assembler}_candidates.tsv")
        cfg = config(assembler, args.bam_file, candidate_sites_file, args.data_name)

        if os.path.exists(cfg.tss_output_file) and os.path.exists(cfg.tes_output_file):
            print(f"✅ Features already extracted for {assembler}. Skipping...")
            continue

        # Extract features
        print(f"⏳ Extracting features for {assembler}...")
        extract_features.main(cfg)
        print(f"✅ Feature extraction complete for {assembler}!")
    
    
    print("🔍 Feature extraction completed for all assemblers.")
    label_cmd = ['python', 'src/label_candidates.py', 
                 '-f', feature_dir,
                 '-o', train_dir,
                 '-r', reference_file,
                 '-m', chrom_mapping,
                 '-d', '50',
                 '-b']
    # Label candidates

    out = subprocess.run(label_cmd)
    # print(out.stdout)
    # print(out.stderr)
    if out.returncode != 0:
        print("❌ Error during candidate labeling!")
        sys.exit(1)
    else:
        print("✅ Candidate labeling complete!")
    
    # Train models
    out_dir = f"out/{args.data_name}"
    os.makedirs(out_dir, exist_ok=True)
    log_dir = f"logs/{args.data_name}"
    os.makedirs(log_dir, exist_ok=True)
    print("⏳ Training models...")
    train_all.main(log_dir, train_dir, out_dir)
    print("✅ Model training complete!")


    # Generate ROC data
    roc_out_dir = generate_roc_data.main(args.data_name)
    plot_out_dir = f"out/{args.data_name}/plots"
    # Plot PR curves
    for assembler in assemblers:
        plot_pr_curves(roc_out_dir, plot_out_dir, assembler)


if __name__ == "__main__":
    main()