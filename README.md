# Telos

Telos scores transcript start sites (TSS), transcript end sites (TES), and assembled
transcripts from a coordinate-sorted BAM and a transcript GTF. It works with any
assembler that can provide a numeric transcript abundance (StringTie, IsoQuant,
Scallop2, …).

```text
# Predict (pretrained or your own models)
BAM + assembly.gtf + models/  →  telos predict  →  scored / ranked transcripts

# Train
BAM + assembly.gtf + reference.gtf + assembly.tmap  →  telos train  →  models/
```

## Install

```bash
git clone https://github.com/Shao-Group/Telos.git
cd Telos
conda env create -f environment.yml
conda activate telos
pip install -e . --no-deps
telos --version   # expect: telos 2.0.0
```

This installs Python dependencies **and** `gffcompare` (needed only to build a `.tmap`
for training, or when using `telos prepare-gtf --make-tmap`).

Requires **Python ≥ 3.10**.

When the Bioconda package is published:

```bash
conda install -c bioconda -c conda-forge telos
```

## What Telos expects

| Command | Required inputs |
| --- | --- |
| `telos predict` | Sorted BAM + index (`.bai`/`.csi`), assembly GTF with transcript `cov "..."`, model directory |
| `telos train` | Same BAM/GTF, plus reference annotation GTF and a gffcompare `.tmap` |
| `telos prepare-gtf` | Assembler GTF; optional IsoQuant TPM table; optional `--make-tmap` + reference GTF |

Model directory (8 files from `telos train` or a pretrained bundle):

```text
models/
  stage1_tss_{rf,xgb}_model.joblib
  stage1_tes_{rf,xgb}_model.joblib
  stage2_model_{rf,xgb}.joblib
  stage2_feature_names_{rf,xgb}.json
```

## Prepare inputs (IsoQuant and tmap)

Telos reads abundance from the transcript attribute `cov "..."`. Assemblers differ:

| Assembler | Typical prep |
| --- | --- |
| **StringTie / Scallop2** | Already have `cov`. For training only, build a tmap. |
| **IsoQuant** | Abundance is often in `*.transcript_model_tpm.tsv`. Inject it into `cov`, then (for training) build a tmap. |

### IsoQuant → Telos-ready GTF

```bash
telos prepare-gtf \
  --gtf OUT.transcript_models.gtf \
  --tpm OUT.transcript_model_tpm.tsv \
  --outdir prepared
# writes prepared/assembly.gtf
```

### Build a tmap for training

```bash
# StringTie / Scallop2 (already have cov)
telos prepare-gtf \
  --gtf stringtie.gtf \
  --ref-gtf gencode.gtf \
  --make-tmap \
  --outdir prepared

# IsoQuant (TPM → cov, then tmap)
telos prepare-gtf \
  --gtf OUT.transcript_models.gtf \
  --tpm OUT.transcript_model_tpm.tsv \
  --ref-gtf gencode.gtf \
  --make-tmap \
  --outdir prepared
# writes prepared/assembly.gtf and prepared/assembly.tmap
```

Use those paths with `telos train` / `telos predict`.

## Quick start: predict

```bash
telos predict \
  --bam sample.sorted.bam \
  --gtf prepared/assembly.gtf \
  --model-dir /path/to/models \
  --backend xgb \
  --outdir predict_out \
  --stage1-workers 8 \
  -v
```

Useful options:

- `--backend {xgb,rf,both}` — which backend writes the scored GTF (default `xgb`). Both
  backends are always written to ranked TSVs.
- `--min-score 0.8` — also write a filtered high-confidence GTF.
- `--config FILE` — use the same Stage I config used during training if you customized it.

### Prediction outputs

| Path | Contents |
| --- | --- |
| `predictions/transcripts.scored.<backend>.gtf` | Input GTF with transcript `cov` set to the Telos score |
| `predictions/transcripts.filtered.<backend>.gtf` | Score-filtered GTF (only with `--min-score`) |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | Ranked `transcript_id`, `pred_prob`, `pred_label` |
| `predictions/sites.scored.tsv` | Site-level RF/XGB probabilities |
| `reports/summary.txt` | Output paths (also printed at completion) |

## Train your own models

```bash
telos train \
  --bam sample.sorted.bam \
  --gtf prepared/assembly.gtf \
  --ref-gtf gencode.gtf \
  --tmap prepared/assembly.tmap \
  --outdir train_out \
  --stage1-workers 8 \
  -v
```

Then predict with `train_out/models/`.

### Training outputs

| Path | Contents |
| --- | --- |
| `models/…` | Complete 8-file model bundle |
| `predictions/sites.scored.tsv` | Site scores on the training assembly |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | Transcript rankings (+ diagnostic labels) |
| `reports/train_metrics.csv` | Stage I / Stage II validation metrics |

## Pretrained models

Model weights are large and are **not** stored in git. Download a modality bundle
(GitHub Release / Zenodo — links TBD) and point `--model-dir` at its `models/` folder.

| Modality | Train sample | Assemblers | Use for |
| --- | --- | --- | --- |
| `ont-cdna` | ENCFF023EXJ | StringTie + IsoQuant | ONT cDNA |
| `ont-drna` | NA12878 | StringTie + IsoQuant | ONT direct RNA |
| `pacbio` | ENCFF450VAU | StringTie + IsoQuant | PacBio Iso-Seq |
| `short-reads` | SRR307903 | StringTie + Scallop2 | Illumina RNA-seq |

```bash
telos predict \
  --bam your.sorted.bam \
  --gtf your.assembly.gtf \
  --model-dir telos-models-ont-cdna/models \
  --outdir predict_out
```

Pick the modality closest to your library type. Prefer training on matched data when
your protocol or species differs substantially.

Paper benchmarks and reproduction live in
[Telos-test](https://github.com/Shao-Group/Telos-test).

## License and citation

Telos is available under the BSD 3-Clause License.

```bibtex
@article{Telos,
  title = {Boosting Transcript Assembly via Delineating Transcript Start and End Sites},
  url = {http://dx.doi.org/10.1101/2025.10.13.682211},
  DOI = {10.1101/2025.10.13.682211},
  publisher = {Cold Spring Harbor Laboratory},
  author = {Khan, Irtesam Mahmud and Zang, Xiaofei Carl and Teng, Ange and Zahin, Tasfia and Shao, Mingfu},
  year = {2025},
  month = oct
}
```
