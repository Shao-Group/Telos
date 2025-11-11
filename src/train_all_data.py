import pandas as pd
from install import install
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

RNASEQ_DIR = "../../tools/rnaseqtools"
OUTPUT_DIR = "train_output"
GENCODE_REF = "data/GRCh38_gencode.gtf"
ENSEMBLE_REF = "data/GRCh38_ensembl.gtf"
PROJECT_CONFIG_DIR = "project_config"
MODEL_CONFIG_DIR = "project_config"
GFFCOMPARE_ENV = "gffcompare"
LOG_DIR = "logs"

def train_data(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file):
    output_dir = os.path.join(output_dir, prefix)
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(LOG_DIR, f"{prefix}.log")
    with open(log_file, "w") as log:
        log.write(f"====== Processing {prefix}... ========\n")
        print(f"Installing {prefix}...")
        config_file_path = install(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file)
        print(f"Installation of {prefix} completed.")

        
        
        # Update config to use enhanced features for training
        print(f"Configuring {prefix} to use enhanced features for training...")
        from config import load_config, save_config
        temp_cfg = load_config(config_file_path)
        temp_cfg.use_embeddings = True
        temp_cfg.embedding_type = 'kmer'
        temp_cfg.embedding_mode = 'hybrid'  # Use k-mer embeddings + basic features
        temp_cfg.is_training = True  # Enable feature selection during training
        save_config(config_file_path)
        print(f"✅ Configuration updated for enhanced features and training mode.")


        print(f"Running extract_features.py for {prefix}...")
        p = subprocess.run(
            ["python", "src/extract_features.py", "--config", config_file_path],
             capture_output=True, text=True
        )
        log.write(p.stdout)
        log.write(p.stderr)
        log.flush()
        if p.returncode != 0:
            print(f"❌ Error in feature extraction for {prefix}: {p.stderr}")
            exit(1)
        else:
            print(f"✅ Feature extraction for {prefix} completed.")
        
        print(f"Running label_candidates.py for {prefix}...")
        p = subprocess.run(
            ["python", "src/label_candidates.py", "--config", config_file_path],
            capture_output=True, text=True, check=True
        )
        log.write(p.stdout)
        log.write(p.stderr)

        print(f"Training model for {prefix}...")
        p = subprocess.run(
            ["python", "src/train_model.py", "--project-config", config_file_path, "--model-config-folder", MODEL_CONFIG_DIR],
            capture_output=True, text=True
        )
        log.write(p.stdout)
        log.write(p.stderr)
        log.flush()
        if p.returncode != 0:
            print(f"❌ Error in model training for {prefix}: {p.stderr}")
            exit(1)
        else:
            print(f"✅ Model training for {prefix} completed.")

        print(f"Generating ROC data for {prefix}...")
        p = subprocess.run(
            ["python", "src/generate_roc_data.py", "--project-config", config_file_path, "--gffcompare-env", GFFCOMPARE_ENV],
            capture_output=True, text=True
        )
        log.write(p.stdout)
        log.write(p.stderr)
        log.flush()
        if p.returncode != 0:
            print(f"❌ Error in ROC data generation for {prefix}: {p.stderr}")
            exit(1)
        else:
            print(f"✅ ROC data generation for {prefix} completed.")

def process_dataset_tools(row_data, rnaseq_dir, output_dir):
    """Process both tools for a single dataset in parallel"""
    prefix, bam_file, gtf_file1, gtf_file2, tmap_file1, tmap_file2, ref = row_data
    
    print(f"Processing {prefix}...")
    prefix1 = prefix + "_stringtie"
    prefix2 = (prefix + "_isoquant") if not prefix.startswith("SRR") else (prefix + "_scallop2")
    
    # Create a list of tasks for this dataset
    tasks = [
        (prefix1, rnaseq_dir, output_dir, bam_file, gtf_file1, ref, tmap_file1),
        (prefix2, rnaseq_dir, output_dir, bam_file, gtf_file2, ref, tmap_file2)
    ]
    
    results = []
    for task in tasks:
        results.append(train_data(*task))
                
    return results

def init_parallel():
    """Parallel version of init() function"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    train_configs = {
        # "prefix": ["cDNA-NA12878","dRNA-NA12878", "pacbio_ENCFF450VAU", "SRR307903"],
        # "prefix" : ["cDNA-NA12878","dRNA-ENCFF155CFF", "pacbio_ENCFF450VAU", "SRR307903", "cDNA-ENCFF023EXJ", "dRNA-NA12878"],
        "prefix" : ["cDNA-ENCFF023EXJ","dRNA-NA12878", "pacbio_ENCFF450VAU", "SRR307903"],
        "bam_file": [
                     "data/nanopore_cDNA_ENCFF023EXJ/ENCFF023EXJ.sorted.bam",
                     "data/nanopore_dRNA_NA12878/NA12878-DirectRNA.sorted.bam",
                     "data/pacbio_ENCFF450VAU/ENCFF450VAU.sorted.bam",
                     "data/SRR307903_hisat/hisat.sorted.bam"
                    ],
                    #  "data/nanopore_cDNA_ENCFF023EXJ/ENCFF023EXJ.sorted.bam",
                    #  "data/nanopore_dRNA_NA12878/NA12878-DirectRNA.sorted.bam"],
        "gtf_file1": [
        # ["data/nanopore_cDNA_NA12878/stringtie_filtered.gtf", 
                      "data/nanopore_cDNA_ENCFF023EXJ/stringtie_filtered.gtf",
                      "data/nanopore_dRNA_NA12878/stringtie_filtered.gtf",
                      "data/pacbio_ENCFF450VAU/stringtie_filtered.gtf",
                      "data/SRR307903_hisat/stringtie_filtered.gtf"],
                    #   "data/nanopore_cDNA_ENCFF023EXJ/stringtie_filtered.gtf",
                    #   "data/nanopore_dRNA_NA12878/stringtie_filtered.gtf"],
        "gtf_file2": ["data/nanopore_cDNA_ENCFF023EXJ/isoquant_filtered.gtf",
                      "data/nanopore_dRNA_NA12878/isoquant_filtered.gtf",
                      "data/pacbio_ENCFF450VAU/isoquant_filtered.gtf",
                      "data/SRR307903_hisat/scallop2_filtered.gtf"],
                    #   "data/nanopore_cDNA_NA12878/isoquant_filtered.gtf", ,
                    #   "data/nanopore_dRNA_NA12878/isoquant_filtered.gtf"],
        "tmap_file1":  ["data/nanopore_cDNA_ENCFF023EXJ/stringtie.stringtie_filtered.gtf.tmap",
                        "data/nanopore_dRNA_NA12878/stringtie.stringtie_filtered.gtf.tmap",
                       "data/pacbio_ENCFF450VAU/stringtie.stringtie_filtered.gtf.tmap",
                       "data/SRR307903_hisat/stringtie.stringtie_filtered.gtf.tmap"],
                    #    "data/nanopore_cDNA_NA12878/stringtie.stringtie_filtered.gtf.tmap",
                    #    "data/nanopore_dRNA_NA12878/stringtie.stringtie_filtered.gtf.tmap"],
        "tmap_file2":  ["data/nanopore_cDNA_ENCFF023EXJ/isoquant.isoquant_filtered.gtf.tmap",   
                        "data/nanopore_dRNA_NA12878/isoquant.isoquant_filtered.gtf.tmap",
                       "data/pacbio_ENCFF450VAU/isoquant.isoquant_filtered.gtf.tmap",
                        "data/SRR307903_hisat/scallop2.scallop2_filtered.gtf.tmap"],
                    #    "data/nanopore_cDNA_NA12878/isoquant.isoquant_filtered.gtf.tmap",
                    #    "data/nanopore_dRNA_NA12878/isoquant.isoquant_filtered.gtf.tmap"],
        "ref_anno_gtf": [GENCODE_REF, GENCODE_REF, GENCODE_REF, ENSEMBLE_REF] #, GENCODE_REF, GENCODE_REF]
    }
    
    train_configs_df = pd.DataFrame(train_configs)
    
    # Prepare data for parallel processing
    dataset_tasks = []
    for index, row in train_configs_df.iterrows():
        row_data = (row["prefix"], row["bam_file"], row["gtf_file1"], 
                   row["gtf_file2"], row["tmap_file1"], row["tmap_file2"], 
                   row["ref_anno_gtf"])
        dataset_tasks.append(row_data)
    
    # dataset_tasks = [dataset_tasks[1]] #[:5]
    # Process datasets in parallel (max 2 concurrent datasets to avoid overwhelming system)
    max_concurrent_datasets = min(4, len(dataset_tasks))
    print(f"🚀 Starting parallel processing with {max_concurrent_datasets} concurrent datasets...")
    
    with ProcessPoolExecutor(max_workers=max_concurrent_datasets) as executor:
        # Submit all dataset processing tasks
        future_to_dataset = {
            executor.submit(process_dataset_tools, task, RNASEQ_DIR, OUTPUT_DIR): task[0]
            for task in dataset_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_dataset):
            dataset_name = future_to_dataset[future]
            try:
                results = future.result()
                print(f"📊 Dataset {dataset_name} processing completed:")
                for result in results:
                    print(f"  {result}")
            except Exception as exc:
                print(f"❌ Dataset {dataset_name} failed with exception: {exc}")

# def init():
#     os.makedirs(OUTPUT_DIR, exist_ok=True)
#     train_configs = {
#         "prefix": ["cDNA-NA12878","dRNA-NA12878", "pacbio_ENCFF450VAU", "SRR307903"],
#         "bam_file": ["data/nanopore_cDNA_NA12878/NA12878-cDNA.sorted.bam",
#                      "data/nanopore_dRNA_NA12878/NA12878-DirectRNA.sorted.bam",
#                      "data/pacbio_ENCFF450VAU/ENCFF450VAU.sorted.bam",
#                      "data/SRR307903_hisat/hisat.sorted.bam"],
#         "gtf_file1": ["data/nanopore_cDNA_NA12878/stringtie.gtf", 
#                       "data/nanopore_dRNA_NA12878/stringtie.gtf",
#                       "data/pacbio_ENCFF450VAU/stringtie.gtf",
#                       "data/SRR307903_hisat/stringtie.gtf"],
#         "gtf_file2": ["data/nanopore_cDNA_NA12878/isoquant.gtf", 
#                       "data/nanopore_dRNA_NA12878/isoquant.gtf",
#                       "data/pacbio_ENCFF450VAU/isoquant.gtf",
#                       "data/SRR307903_hisat/scallop2.gtf"],
#         "tmap_file1":  ["data/nanopore_cDNA_NA12878/stringtie.stringtie.gtf.tmap",
#                         "data/nanopore_dRNA_NA12878/stringtie.stringtie.gtf.tmap",
#                        "data/pacbio_ENCFF450VAU/stringtie.stringtie.gtf.tmap",
#                        "data/SRR307903_hisat/stringtie.stringtie.gtf.tmap"],
#         "tmap_file2":  ["data/nanopore_cDNA_NA12878/isoquant.isoquant.gtf.tmap",
#                         "data/nanopore_dRNA_NA12878/isoquant.isoquant.gtf.tmap",
#                        "data/pacbio_ENCFF450VAU/isoquant.isoquant.gtf.tmap",
#                        "data/SRR307903_hisat/scallop2.scallop2.gtf.tmap"],
#         "ref_anno_gtf": [GENCODE_REF, GENCODE_REF, GENCODE_REF, ENSEMBLE_REF]

#     }
#     train_configs_df = pd.DataFrame(train_configs)
#     for index, row in train_configs_df.iterrows():
#         prefix = row["prefix"]
#         bam_file = row["bam_file"]
#         gtf_file1 = row["gtf_file1"]
#         gtf_file2 = row["gtf_file2"]
#         tmap_file1 = row["tmap_file1"]
#         tmap_file2 = row["tmap_file2"]
#         ref = row["ref_anno_gtf"]

#         print(f"Processing {prefix}...")
#         prefix1 = (prefix + "_stringtie") 
#         prefix2 = (prefix + "_isoquant") if not prefix.startswith("SRR") else (prefix + "_scallop2")
#         # def install(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file):
#         train_data(prefix1, RNASEQ_DIR, OUTPUT_DIR, bam_file, gtf_file1, ref, tmap_file1)
#         train_data(prefix2, RNASEQ_DIR, OUTPUT_DIR, bam_file, gtf_file2, ref, tmap_file2)


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    
    # Use the parallel version by default
    # Change to init() if you want sequential processing
    init_parallel()
    

if __name__ == "__main__":
    main()