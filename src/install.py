import argparse
import os
from config import create_config, reset_config, save_config


def install(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file):
    """
    Install the package.
    """
    print("Initializing Telos ...")
    if not os.path.exists(rnaseq_dir):
        raise ValueError(f"RNAseq directory {rnaseq_dir} does not exist.")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    else:
        print(f"Warning: Output directory {output_dir} already exists.")
        # raise ValueError(f"Output directory {output_dir} already exists. Please choose a different name.")
    
    if not os.path.exists(bam_file):
        raise ValueError(f"BAM file {bam_file} does not exist.")
    if not os.path.exists(gtf_file):
        raise ValueError(f"GTF file {gtf_file} does not exist.")
    
    reset_config()  # Reset the configuration to ensure a clean state
    create_config(
        bam_file=bam_file,
        prefix=prefix,
        output_dir=output_dir,
        rnaseqtools_dir=rnaseq_dir,
        ref_anno_gtf=ref_anno_gtf,
        gtf_file=gtf_file,
        tmap_file=tmap_file
    )
    print("Configuration created.")
    print(f"Saving configuration to 'project_config/{prefix}_config.pkl'...")
    config_file = f"project_config/{prefix}_config.pkl"
    save_config(config_file)
    print(f"Configuration saved to {config_file}.")

    return config_file
    # Add your installation logic here


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Install RNAseq package.")
    parser.add_argument(
        "--dir-rnaseq",
        type=str,
        required=True,
        help="Directory containing the RNAseq package."
    )
    parser.add_argument(
        "--dir-output",
        type=str,
        required=True,
        help="Directory to store output files."
    )
    parser.add_argument(
        "--file-bam",
        type=str,
        required=True,
        help="Path to the input BAM file."
    )
    parser.add_argument(
        "--file-gtf",
        type=str,
        required=True,
        help="Path to the input GTF file."
    )
    parser.add_argument(
        "--prefix",
        type=str,
        required=True,
        help="Prefix for output files."
    )
    parser.add_argument(
        "--ref-anno-gtf",
        type=str,
        required=True,
        help="Path to the reference annotation GTF file."
    )
    parser.add_argument(
        "--tmap-file",
        type=str,
        required=True,
        help="Path to the tmap file."
    )

    args = parser.parse_args()
    install(
        prefix=args.prefix,
        rnaseq_dir=args.dir_rnaseq,
        output_dir=args.dir_output,
        bam_file=args.file_bam,
        gtf_file=args.file_gtf,
        ref_anno_gtf=args.ref_anno_gtf,
        tmap_file=args.tmap_file
    )