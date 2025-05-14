# Find-TSS-TES

This project identifies **Transcript Start Sites (TSS)** and **Transcript End Sites (TES)** from Oxford Nanopore long-read RNA sequencing data, using features extracted from BAM alignments. It is **annotation-agnostic** for training and supports multiple tools and models for benchmarking.

## 🔍 Overview

* Input: BAM files aligned with Minimap2 + candidate sites from StringTie or IsoQuant
* Output: Trained classifiers for TSS and TES, with evaluation metrics and feature importance
* Features extracted: coverage shifts, soft-clipping, entropy, splice junction distance, etc.
* Models supported: XGBoost, LightGBM, RandomForest
* Evaluation: Precision/Recall, F1, AUPR, Accuracy, Confusion Matrix

## 📂 Directory Structure

```
├── features/               # Raw feature CSVs per tool and site type
├── data_train/       # Same features but labeled using RefSeq
├── reports/                # Model evaluation metrics and plots
├── models/                 # Trained model artifacts
├── configs/                # YAML config files for TSS/TES models
├── logs/                   # CLI logs and benchmark summaries
├── src/                    # All scripts and utilities
├── out/                    # All output files 
```

## ⚙️ Scripts & How to Run

### 1. Extract Features from BAM

Extract features for candidate sites from a BAM file.

```bash
python src/extract_features.py --method stringtie (or isoquant) 
```

### 2. Label Features using Reference Annotation

Assign labels (1 = true site, 0 = false) using a distance threshold to reference TSS/TES.

For batch labeling, 
```bash
python src/label_candidates.py \
  --reference data/refSeq.tsstes \
  --distance 50 \
  --mapping data/GRCh38_RefSeq2UCSC.txt \
  --batch 
```

### 3. Train Models for a Given Tool

Trains models for both TSS and TES using a specific tool's candidates.

```bash
python src/train_all.py 
```


### 4. Summarize Metrics

All the stats, summarization and plot generation available at src/data-analysis.ipynb



### Data Names

  - nanopore_cDNA_NA12878  
  - pacbio_ENCFF450VAU  
  - nanopore_dRNA_NA12878   
  - pacbio_ENCFF694DIE  
  - SRR307903_hisat
---

## ✍️ Author

Developed by [irtesampsu](https://github.com/irtesampsu) as part of course project for CSE 566 at Penn State.

---

For help or issues, open an issue on GitHub or contact the author.
