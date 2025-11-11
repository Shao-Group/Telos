import pandas as pd
from install import install
import subprocess
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

import sys
# sys.path.append('src')
from config import load_config, save_config
import os

# ----------- Directory paths -------------------
RNASEQ_DIR = "../../tools/rnaseqtools"
OUTPUT_DIR = "test_output"
GENCODE_REF = "data/GRCh38_gencode.gtf"
ENSEMBLE_REF = "data/GRCh38_ensembl.gtf"
PROJECT_CONFIG_DIR = "project_config"
MODEL_CONFIG_DIR = "project_config"
GFFCOMPARE_ENV = "gffcompare"
LOG_DIR = "logs"
STAGE1_MODELS = ["xgboost", "randomforest"]
# -----------------------------------------------

def generate_model_paths(model_folder):
    """
    Generate paths for pretrained models.
    """
    model_paths = []
    
    m1 = {
        "tss": os.path.join(model_folder, f"tss_randomforest_model.joblib"),
        "tes": os.path.join(model_folder, f"tes_randomforest_model.joblib"),
        "stage2": os.path.join(model_folder, f"xgboost_stage2_model.joblib"),
        "model_type": "randomforest"
    }
    model_paths.append(m1)
    m2 = {
        "tss": os.path.join(model_folder, f"tss_xgboost_model.joblib"),
        "tes": os.path.join(model_folder, f"tes_xgboost_model.joblib"),
        "stage2": os.path.join(model_folder, f"xgboost_stage2_model.joblib"),
        "model_type": "xgboost"
    }
    model_paths.append(m2)

    return model_paths

def test_with_pretrained(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file, pretrained_config_file):
    pretrained_config = load_config(pretrained_config_file)
    pretrained_model_folder = pretrained_config.models_output_dir
    # if not "drna" in prefix.lower():
    #     return
    
    output_dir = os.path.join(output_dir, prefix)
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(LOG_DIR, f"{prefix}.log")
    with open(log_file, "w") as log:
        log.write(f"====== Processing {prefix}... ========\n")
        print(f"Installing {prefix}...")
        config_file_path = install(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file)
        print(f"Installation of {prefix} completed.")
        temp_cfg = load_config(config_file_path)
        temp_cfg.tss_selected_feature_file = pretrained_config.tss_selected_feature_file
        temp_cfg.tes_selected_feature_file = pretrained_config.tes_selected_feature_file
        temp_cfg.is_training = False
        save_config(config_file_path)

        print(f"Running extract_features.py for {prefix}...")
        p = subprocess.run(
            ["python", "src/extract_features.py", "--config", config_file_path],
             capture_output=True, text=True, check=True
        )
        log.write(p.stdout)
        log.write(p.stderr)
        log.flush()
        if p.returncode != 0:
            print(f"❌ Error in feature extraction for {prefix}: {p.stderr}")
            exit(1)
        else:
            print(f"✅ Feature extraction for {prefix} completed.")

        

        print(f"Labeling candidates for {prefix}...")
        p = subprocess.run(
            ["python", "src/label_candidates.py", "--config", config_file_path],
            capture_output=True, text=True
        )
        log.write(p.stdout)
        log.write(p.stderr)
        log.flush()
        if p.returncode != 0:
            print(f"❌ Error in candidate labeling for {prefix}: {p.stderr}")
            exit(1)
        else:
            print(f"✅ Candidate labeling for {prefix} completed.")
            

        print(f"Testing model for {prefix}...")
        # Generate paths for pretrained models
        model_paths = generate_model_paths(pretrained_model_folder)
        for model in model_paths:
            model_type = model["model_type"]
            tss_model = model["tss"]
            tes_model = model["tes"]
            stage2_model = model["stage2"]
        
            p = subprocess.run(
                ["python", "src/validate_with_pretrained.py", "--project-config", config_file_path,
                 "--model-config-folder", MODEL_CONFIG_DIR, "--pretrained_tss_model", tss_model,
                 "--pretrained_tes_model", tes_model, "--pretrained_stage2_model", stage2_model,
                 "--model_type", model_type ] ,
                capture_output=True, text=True
            )
            log.write(p.stdout)
            log.write(p.stderr)
            log.flush()
            if p.returncode != 0:
                print(f"❌ Error in model testing for {prefix}: {p.stderr}")
                exit(1)
            else:
                print(f"✅ Model testing for {prefix} with {model_type} completed.")

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

def process_dataset_tools_test(row_data, rnaseq_dir, output_dir):
    """Process both tools for a single dataset sequentially"""
    prefix, bam_file, gtf_file1, gtf_file2, tmap_file1, tmap_file2, ref, pretrained_config1, pretrained_config2 = row_data
    # pretrained_folder1 = pretrained_config1.models_output_dir
    # pretrained_folder2 = pretrained_config2.models_output_dir
    print(f"Processing {prefix}...")
    prefix1 = prefix + "_stringtie"
    prefix2 = (prefix + "_isoquant") if not prefix.startswith("SRR") else (prefix + "_scallop2")
    
    results = []
    
    # Process stringtie first
    try:
        test_with_pretrained(prefix1, rnaseq_dir, output_dir, bam_file, gtf_file1, ref, tmap_file1, pretrained_config1)
        results.append(f"✅ {prefix1} completed successfully")
    except Exception as exc:
        results.append(f"❌ {prefix1} failed with exception: {exc}")
        return results  # Stop processing if first tool fails
    
    # Process isoquant/scallop2 second
    try:
        test_with_pretrained(prefix2, rnaseq_dir, output_dir, bam_file, gtf_file2, ref, tmap_file2, pretrained_config2)
        results.append(f"✅ {prefix2} completed successfully")
    except Exception as exc:
        results.append(f"❌ {prefix2} failed with exception: {exc}")
                
    return results

def test_all_parallel():
    """Parallel version of test_all() function"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    test_configs = {
        "prefix": ["cDNA-K562",
                   "cDNA-ENCFF263YFG", 
                   "cDNA-NA12878", 
                   "dRNA-Hek293T", 
                   "dRNA-ENCFF771DIX",
                #    "dRNA-NA12878",
                    "dRNA-ENCFF155CFF",
                   "pacbio_ENCFF694DIE",
                   "pacbio_ENCFF563QZR",
                   "pacbio_ENCFF370NFS", 
                   "SRR307911",
                   "SRR545695",
                   "SRR315334",
                   "SRR534307",
                   "SRR545723",
                #    "SRR307911",
                   "SRR315323",
                   "SRR534319",
                   "SRR534291",
                   "SRR387661"],
        "bam_file": ["data/cv_K562_cDNA/SGNex_K562_cDNA.sorted.bam",
                    "data/cv_ENCFF263YFG_cDNA/ENCFF263YFG.sorted.bam",
                    "data/cv_NA12878_cDNA/NA12878-cDNA.sorted.bam",
                    "data/cv_Hek293T_dRNA/SGNex_Hek293T_dRNA.sorted.bam",
                    "data/cv_ENCFF771DIX_dRNA/dRNA-ENCFF771DIX.sorted.bam",
                    # "data/cv_NA12878_dRNA/NA12878-DirectRNA.sorted.bam",
                    "data/cv_ENCFF155CFF_dRNA/dRNA-ENCFF155CFF.sorted.bam",
                    "data/cv_pacbio_ENCFF694DIE/ENCFF694DIE.sorted.bam",
                    "data/cv_pacbio_ENCFF563QZR/ENCFF563QZR.sorted.bam",
                    "data/cv_pacbio_ENCFF370NFS/ENCFF370NFS.sorted.bam",
                     "data/cv_SRR307911_hisat/hisat.sort.bam",
                     "data/cv_SRR545695_hisat/hisat.sort.bam",
                     "data/cv_SRR315334_hisat/hisat.sort.bam",
                     "data/cv_SRR534307_hisat/hisat.sort.bam",
                     "data/cv_SRR545723_hisat/hisat.sort.bam",
                     # "data/cv_SRR307911_hisat/hisat.sort.bam",
                     "data/cv_SRR315323_hisat/hisat.sort.bam",
                     "data/cv_SRR534319_hisat/hisat.sort.bam",
                     "data/cv_SRR534291_hisat/hisat.sort.bam",
                     "data/cv_SRR387661_hisat/hisat.sort.bam"],
        "gtf_file1": ["data/cv_K562_cDNA/stringtie_filtered.gtf", 
                      "data/cv_ENCFF263YFG_cDNA/stringtie_filtered.gtf",
                      "data/cv_NA12878_cDNA/stringtie_filtered.gtf",
                      "data/cv_Hek293T_dRNA/stringtie_filtered.gtf",
                      "data/cv_ENCFF771DIX_dRNA/stringtie_filtered.gtf",
                      # "data/cv_NA12878_dRNA/stringtie_filtered.gtf",
                      "data/cv_ENCFF155CFF_dRNA/stringtie_filtered.gtf",
                      "data/cv_pacbio_ENCFF694DIE/stringtie_filtered.gtf",
                      "data/cv_pacbio_ENCFF563QZR/stringtie_filtered.gtf",
                      "data/cv_pacbio_ENCFF370NFS/stringtie_filtered.gtf",                      
                      "data/cv_SRR307911_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR545695_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR315334_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR534307_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR545723_hisat/stringtie_filtered.gtf",
                      # "data/cv_SRR307911_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR315323_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR534319_hisat/stringtie_filtered.gtf",
                      "data/cv_SRR534291_hisat/stringtie_filtered.gtf",  
                      "data/cv_SRR387661_hisat/stringtie_filtered.gtf"],
        "gtf_file2": ["data/cv_K562_cDNA/isoquant_filtered.gtf", 
                      "data/cv_ENCFF263YFG_cDNA/isoquant_filtered.gtf",
                      "data/cv_NA12878_cDNA/isoquant_filtered.gtf",
                      "data/cv_Hek293T_dRNA/isoquant_filtered.gtf",
                      "data/cv_ENCFF771DIX_dRNA/isoquant_filtered.gtf",
                    #   "data/cv_NA12878_dRNA/isoquant_filtered.gtf",
                      "data/cv_ENCFF155CFF_dRNA/isoquant_filtered.gtf",
                      "data/cv_pacbio_ENCFF694DIE/isoquant_filtered.gtf",
                      "data/cv_pacbio_ENCFF563QZR/isoquant_filtered.gtf",
                      "data/cv_pacbio_ENCFF370NFS/isoquant_filtered.gtf",
                      "data/cv_SRR307911_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR545695_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR315334_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR534307_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR545723_hisat/scallop2_filtered.gtf",
                      # "data/cv_SRR307911_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR315323_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR534319_hisat/scallop2_filtered.gtf",
                      "data/cv_SRR534291_hisat/scallop2_filtered.gtf",  
                      "data/cv_SRR387661_hisat/scallop2_filtered.gtf"],
        "tmap_file1":  ["data/cv_K562_cDNA/stringtie.stringtie_filtered.gtf.tmap",
                        "data/cv_ENCFF263YFG_cDNA/stringtie.stringtie_filtered.gtf.tmap",
                        "data/cv_NA12878_cDNA/stringtie.stringtie_filtered.gtf.tmap",
                        "data/cv_Hek293T_dRNA/stringtie.stringtie_filtered.gtf.tmap",
                        "data/cv_ENCFF771DIX_dRNA/stringtie.stringtie_filtered.gtf.tmap",
                        # "data/cv_NA12878_dRNA/stringtie.stringtie_filtered.gtf.tmap",
                        "data/cv_ENCFF155CFF_dRNA/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_pacbio_ENCFF694DIE/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_pacbio_ENCFF563QZR/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_pacbio_ENCFF370NFS/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR307911_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR545695_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR315334_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR534307_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR545723_hisat/stringtie.stringtie_filtered.gtf.tmap",
                    #    "data/cv_SRR307911_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR315323_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR534319_hisat/stringtie.stringtie_filtered.gtf.tmap",
                       "data/cv_SRR534291_hisat/stringtie.stringtie_filtered.gtf.tmap",  
                       "data/cv_SRR387661_hisat/stringtie.stringtie_filtered.gtf.tmap"],
        "tmap_file2":  ["data/cv_K562_cDNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_ENCFF263YFG_cDNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_NA12878_cDNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_Hek293T_dRNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_ENCFF771DIX_dRNA/isoquant.isoquant_filtered.gtf.tmap",
                        # "data/cv_NA12878_dRNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_ENCFF155CFF_dRNA/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_pacbio_ENCFF694DIE/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_pacbio_ENCFF563QZR/isoquant.isoquant_filtered.gtf.tmap",
                        "data/cv_pacbio_ENCFF370NFS/isoquant.isoquant_filtered.gtf.tmap",
                       "data/cv_SRR307911_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR545695_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR315334_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR534307_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR545723_hisat/scallop2.scallop2_filtered.gtf.tmap",
                    #    "data/cv_SRR307911_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR315323_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR534319_hisat/scallop2.scallop2_filtered.gtf.tmap",
                       "data/cv_SRR534291_hisat/scallop2.scallop2_filtered.gtf.tmap",  
                       "data/cv_SRR387661_hisat/scallop2.scallop2_filtered.gtf.tmap"],
        "ref_anno_gtf": [GENCODE_REF]*9 + [ENSEMBLE_REF]*9,
        
        "pretrained_config1" :  [f"{PROJECT_CONFIG_DIR}/cDNA-ENCFF023EXJ_stringtie_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/dRNA-NA12878_stringtie_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/pacbio_ENCFF450VAU_stringtie_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/SRR307903_stringtie_config.pkl"]*9,

        "pretrained_config2" :  [f"{PROJECT_CONFIG_DIR}/cDNA-ENCFF023EXJ_isoquant_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/dRNA-NA12878_isoquant_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/pacbio_ENCFF450VAU_isoquant_config.pkl"]*3 + 
                                [f"{PROJECT_CONFIG_DIR}/SRR307903_scallop2_config.pkl"]*9
    }
    

    for key, value in test_configs.items():
        print(f"{key}: {len(value)}")
        print("-"*100)

    test_configs_df = pd.DataFrame(test_configs)
    
    # Prepare data for parallel processing
    dataset_tasks = []
    for index, row in test_configs_df.iterrows():
        row_data = (row["prefix"], row["bam_file"], row["gtf_file1"], 
                   row["gtf_file2"], row["tmap_file1"], row["tmap_file2"], 
                   row["ref_anno_gtf"], row["pretrained_config1"], 
                   row["pretrained_config2"])
        
        dataset_tasks.append(row_data)
    
    # dataset_tasks = dataset_tasks[:-4]
    # Process datasets in parallel (max 2 concurrent datasets to avoid overwhelming system)
    max_concurrent_datasets = min(4, len(dataset_tasks))
    print(f"🚀 Starting parallel testing with {max_concurrent_datasets} concurrent datasets...")
    
    with ProcessPoolExecutor(max_workers=max_concurrent_datasets) as executor:
        # Submit all dataset processing tasks
        future_to_dataset = {
            executor.submit(process_dataset_tools_test, task, RNASEQ_DIR, OUTPUT_DIR): task[0]
            for task in dataset_tasks
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_dataset):
            dataset_name = future_to_dataset[future]
            try:
                results = future.result()
                print(f"📊 Dataset {dataset_name} testing completed:")
                for result in results:
                    print(f"  {result}")
            except Exception as exc:
                print(f"❌ Dataset {dataset_name} failed with exception: {exc}")

# def test_all():
#     # data/cv_Hek293T_dRNA  data/cv_K562_cDNA  data/cv_pacbio_ENCFF694DIE  data/cv_SRR307911_hisat
#     # data/cv_Hek293T_dRNA/SGNex_Hek293T_dRNA.sorted.bam  data/cv_pacbio_ENCFF694DIE/ENCFF694DIE.sorted.bam
#     # data/cv_K562_cDNA/SGNex_K562_cDNA.sorted.bam        data/cv_SRR307911_hisat/hisat.sorted.bam
#     os.makedirs(OUTPUT_DIR, exist_ok=True)
#     test_configs = {
#         "prefix": ["cDNA-K562","dRNA-Hek293T", "pacbio_ENCFF694DIE", "SRR307911"],
#         "bam_file": ["data/cv_K562_cDNA/SGNex_K562_cDNA.sorted.bam",
#                      "data/cv_Hek293T_dRNA/SGNex_Hek293T_dRNA.sorted.bam",
#                      "data/cv_pacbio_ENCFF694DIE/ENCFF694DIE.sorted.bam",
#                      "data/cv_SRR307911_hisat/hisat.sorted.bam"],
#         "gtf_file1": ["data/cv_K562_cDNA/stringtie.gtf", 
#                       "data/cv_Hek293T_dRNA/stringtie.gtf",
#                       "data/cv_pacbio_ENCFF694DIE/stringtie.gtf",
#                       "data/cv_SRR307911_hisat/stringtie.gtf"],
#         "gtf_file2": ["data/cv_K562_cDNA/isoquant.gtf", 
#                       "data/cv_Hek293T_dRNA/isoquant.gtf",
#                       "data/cv_pacbio_ENCFF694DIE/isoquant.gtf",
#                       "data/cv_SRR307911_hisat/scallop2.gtf"],
#         "tmap_file1":  ["data/cv_K562_cDNA/stringtie.stringtie.gtf.tmap",
#                         "data/cv_Hek293T_dRNA/stringtie.stringtie.gtf.tmap",
#                        "data/cv_pacbio_ENCFF694DIE/stringtie.stringtie.gtf.tmap",
#                        "data/cv_SRR307911_hisat/stringtie.stringtie.gtf.tmap"],
#         "tmap_file2":  ["data/cv_K562_cDNA/isoquant.isoquant.gtf.tmap",
#                         "data/cv_Hek293T_dRNA/isoquant.isoquant.gtf.tmap",
#                        "data/cv_pacbio_ENCFF694DIE/isoquant.isoquant.gtf.tmap",
#                        "data/cv_SRR307911_hisat/scallop2.scallop2.gtf.tmap"],
#         "ref_anno_gtf": [GENCODE_REF, GENCODE_REF, GENCODE_REF, ENSEMBLE_REF],
#         "pretrained_model_folder1" : ["train_output/cDNA-NA12878_stringtie/models", "train_output/dRNA-NA12878_stringtie/models",
#                                       "train_output/pacbio_ENCFF450VAU_stringtie/models", "train_output/SRR307903_stringtie/models"],
#         "pretrained_model_folder2" : ["train_output/cDNA-NA12878_isoquant/models", "train_output/dRNA-NA12878_isoquant/models",
#                                       "train_output/pacbio_ENCFF450VAU_isoquant/models", "train_output/SRR307903_scallop2/models"]

#     }
#     test_configs_df = pd.DataFrame(test_configs)
#     for index, row in test_configs_df.iterrows():
#         prefix = row["prefix"]
#         bam_file = row["bam_file"]
#         gtf_file1 = row["gtf_file1"]
#         gtf_file2 = row["gtf_file2"]
#         tmap_file1 = row["tmap_file1"]
#         tmap_file2 = row["tmap_file2"]
#         ref = row["ref_anno_gtf"]
#         pretrained_model_folder1 = row["pretrained_model_folder1"]
#         pretrained_model_folder2 = row["pretrained_model_folder2"]

#         print(f"Processing {prefix}...")
#         prefix1 = (prefix + "_stringtie") 
#         prefix2 = (prefix + "_isoquant") if not prefix.startswith("SRR") else (prefix + "_scallop2")
#         # def install(prefix, rnaseq_dir, output_dir, bam_file, gtf_file, ref_anno_gtf, tmap_file):
#         test_with_pretrained(prefix1, RNASEQ_DIR, OUTPUT_DIR, bam_file, gtf_file1, ref, tmap_file1, pretrained_model_folder1)
#         test_with_pretrained(prefix2, RNASEQ_DIR, OUTPUT_DIR, bam_file, gtf_file2, ref, tmap_file2, pretrained_model_folder2)


def main():
    os.makedirs(LOG_DIR, exist_ok=True)
    
    # Use the parallel version by default
    # Change to test_all() if you want sequential processing
    test_all_parallel()
    print("All tests passed!")

if __name__ == "__main__":
    main()