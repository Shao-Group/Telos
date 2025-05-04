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
```

## ⚙️ Scripts & How to Run

### 1. Extract Features from BAM

Extract features for candidate sites from a BAM file.

```bash
python src/extract_features.py \
  --bam path/to/input.bam \
  --candidates path/to/stringtie_tss_candidates.tsv \
  --output features/stringtie_tss.csv \
  --site_type tss \
  --window_size 100
```

### 2. Label Features using Reference Annotation

Assign labels (1 = true site, 0 = false) using a distance threshold to reference TSS/TES.

```bash
python src/label_candidates.py \
  --candidates features/stringtie_tss.csv \
  --reference refseq_tss.tsv \
  --output features_labeled/stringtie_tss_labeled.csv \
  --site_type tss \
  --distance 50 \
  --mapping chrom_map.tsv
```

### 3. Train Models for a Given Tool

Trains models for both TSS and TES using a specific tool's candidates.

```bash
python src/train_all.py --tool stringtie --model_type xgboost
```

To train across **all tools**, models, and site types:

```bash
python src/train_all.py
```

### 4. Summarize Metrics

Combine all metric summaries from reports/ into one Markdown file.

```bash
python src/generate_metrics_summary.py
```

Output will be saved to `reports/metrics_summary.md`

## 🔬 Feature Engineering

Each candidate site is described by:

* Local read pileup: start/end density, soft-clip stats, full-length counts
* Mapping stats: MAPQ mean/std, strand bias
* Regional structure: coverage change, distance to splice junction
* Entropy: read start/end and soft-clip base entropy

New features like polyA tails and motif scans can be added in future work.

## 📈 Evaluation Outputs

* `metrics_summary.txt`: All key metrics (F1, AUPR, Accuracy, Confusion matrix)
* `pr_curve.png`: Precision-Recall curve per model
* `feature_importance.png`: Top features per model (Random Forest only)

---

## ✍️ Author

Developed by [irtesampsu](https://github.com/irtesampsu) as part of transcript boundary prediction research at Penn State.

---

For help or issues, open an issue on GitHub or contact the author.
