# Telos

Telos is an **annotation-agnostic** tool that improves transcript assembly accuracy by delineating **Transcript Start Sites (TSS)** and **Transcript End Sites (TES)** from RNA sequencing data, using features extracted from BAM alignments. It can be used with any assembler as long as it produces a gtf file annotated with coverage information 

## 🔍 Overview

* Input: BAM files aligned with Minimap2 + assembled transcripts (gtf file) + reference annotation + gffcompare tmap on baseline gtf file.
* Output: Trained classifiers for TSS and TES, and **Transcript Scoring** with evaluation metrics and feature importance
* Features extracted: read density, coverage shifts, soft-clipping, entropy, splice junction distance, etc.
* Models supported: XGBoost, RandomForest
* Evaluation: Precision/Recall, F1, AUPR, Accuracy, Confusion Matrix

## 📂 Directory Structure

```
├── project_config/         # YAML and pkl config files for TSS/TES models and project configurations
├── src/                    # All scripts and utilities
```
Output directory structure:

```
├── data                 # Contains candidate sites, coverage file, and validation baseline
├── features             # Extracted features from the candidate sites
├── models               # Saved models after training
├── predictions          # Stage 1 and 2 prediction on the validation dataset
├── updated_cov          # GTF files after updating the coverage based on predictions
├── reports
  ├── feature_importance # Feature importance from TSS/TES randomforest model
  ├── gffcompare         # GFFCompare results folder
  ├── pr_data            # Precision Recall curve data for stage 1 models
  ├── transcript_pr_data # Transcript level precision recall data
  
```

## ⚙️ Scripts & How to Run
All scripts are run from the **root** folder of the project where this `README.md` is located.

### Step 0: Prerequisites
All the required packages are available at `environment.yml` file. We recommend using anaconda creating a virtual environment using anaconda from the yml file. This can be done using the following command:
```
conda env create -f environment.yml
```
After installing required python packages, [RNASeqtools](https://github.com/Shao-Group/rnaseqtools) needs to be installed following the directive in the source repository. Next, install GFFCompare using bioconda in a separate conda environment in order to run `gffcompare` easily. This is required for generating results comparing Telos with the baseline. You are now ready to run Telos. 
Now, you need the following input files:
  - BAM File: an aligned BAM file of RNA-seq data
  - GTF FILE: Assembled transcripts in a gtf file
  - Reference annotation gtf


### 1. Project setup
Run GFFCompare on the baseline assembled gtf file. You will need to pass the `.tmap` file as input to `install.py`.

Setup the project for analysis.
  - --dir-rnaseq : RNASeq-Tools home directory
  - --prefix: a prefix that will be concatenated to some output files
  - --dir-output: Output directory (should be empty)
  - --file-bam: Path to the BAM file containing aligned reads
  - --file-gtf: Path to the GTF file of assembled transcripts
  - --ref-anno-gtf: Path to the Reference annotation GTF file
  - --tmap-file: *.tmap* file obtained after running GFFCompare on the baseline GTF and reference annotation.

```bash
python src/install.py --dir-rnaseq DIR_RNASEQ --dir-output DIR_OUTPUT --file-bam FILE_BAM --file-gtf FILE_GTF --prefix PREFIX --ref-anno-gtf REF_ANNO_GTF --tmap-file TMAP_FILE
```

The `config.pkl` will be inside `project_config/` directory. You need to pass this filepath in later stages of the analysis. 

### 2. Extract Features

Now, extract features using the bam file. By default, the script does parallel processing with number of parallel processes is set to number of available cpu cores (maximum 8). You can set the number of processes manually by passing the argument `--n-processes=$INTEGER` or you can turn off parallel processing completely by passing `--no-parallel`.

```bash
python src/extract_features.py --config CONFIG_PATH
```

### 3. Label Features using Reference Annotation

Assign labels (1 = true site, 0 = false) using a distance threshold to reference TSS/TES.

For labeling candidate boundaries, 
```bash
python src/label_candidates.py \
  -- config CONFIG_PATH
  --distance 50 
```

### 4. Train Models 

Trains models for both Stage 1 and 2 using a candidate features.

```bash
python src/train_all.py --project-config CONFIG_PATH --model-config-folder project_config
```
Model config folder should contain configuration for the stage 1 models. Example can be found in `project_config` folder.

### 5. Validate only using pretrained model

```bash
python src/validate_with_pretrained.py --project-config PROJECT_CONFIG --model-config-folder MODEL_CONFIG_FOLDER  --pretrained_tss_model PRETRAINED_TSS_MODEL_PATH --pretrained_tes_model PRETRAINED_TES_MODEL_PATH --pretrained_stage2_model PRETRAINED_STAGE2_MODEL_PATH --model_type MODEL_TYPE
```

### 5. Generate comparison data with GTFCuff

To generate results for comparison, run:

```bash
python src/generate_roc_data.py --project-config PROJECT_CONFIG --gffcompare-env GFFCOMPARE_ENV
```


---
## ✍️ Author

Developed by [Shao Group](https://github.com/Shao-Group).
---

For help or issues, open an issue on GitHub or contact the author. 
