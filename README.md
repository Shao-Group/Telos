# Telos

Telos scores transcript start sites (TSS), transcript end sites (TES), and assembled
transcripts from a coordinate-sorted BAM and a transcript GTF. It works with any assembler
that writes numeric transcript-level `cov` attributes (StringTie, IsoQuant, Scallop2, …).

```text
BAM + assembly GTF  →  telos predict  →  scored / ranked transcripts
BAM + GTF + ref + tmap  →  telos train  →  reusable model bundle
```

## Install

**Recommended (from source, Telos 2.0):**

```bash
git clone https://github.com/Shao-Group/Telos.git
cd Telos
conda env create -f environment.yml
conda activate telos
pip install -e . --no-deps
telos --version   # expect: telos 2.0.0
```

Requires **Python ≥ 3.10**. Dependencies are listed in `environment.yml` / `pyproject.toml`
(joblib, lightgbm, numpy, pandas, pysam, pyyaml, scikit-learn, xgboost).

When the Bioconda recipe is published for v2:

```bash
conda install -c bioconda -c conda-forge telos
```

Training and prediction are pure Python and do **not** need gffcompare at runtime. (gffcompare is only needed upstream to produce a `.tmap` for training, and for
optional AUPR evaluation.)

## What you need


| Mode        | Inputs                                                                            |
| ----------- | --------------------------------------------------------------------------------- |
| **Predict** | Sorted BAM + index (`.bai`), assembly GTF with transcript `cov`, model directory  |
| **Train**   | Same BAM/GTF, plus reference annotation GTF and the assembly’s gffcompare `.tmap` |


Model directory (8 files):

```text
models/
  stage1_tss_{rf,xgb}_model.joblib
  stage1_tes_{rf,xgb}_model.joblib
  stage2_model_{rf,xgb}.joblib
  stage2_feature_names_{rf,xgb}.json
```



## Quick start: predict

```bash
telos predict \
  --bam sample.sorted.bam \
  --gtf assembly.gtf \
  --model-dir /path/to/models \
  --backend xgb \
  --outdir predict_out \
  --stage1-workers 8 \
  -v
```

Useful options:

- `--backend {xgb,rf,both}` — which backend writes the scored GTF (default `xgb`). Both
backends are always scored into ranked TSVs.
- `--min-score 0.8` — also write a filtered high-confidence GTF i.e., transcripts passing the minimum score predicted by the model.
- `--config FILE` — use the same Stage I config that was used for training if you customized it.



### Prediction outputs


| Path                                             | Contents                                               |
| ------------------------------------------------ | ------------------------------------------------------ |
| `predictions/transcripts.scored.<backend>.gtf`   | Input GTF with transcript `cov` set to the Telos score |
| `predictions/transcripts.filtered.<backend>.gtf` | Score-filtered GTF (only with `--min-score`)           |
| `predictions/transcripts.ranked.{rf,xgb}.tsv`    | Ranked `transcript_id`, `pred_prob`, `pred_label`      |
| `predictions/sites.scored.tsv`                   | Site-level RF/XGB probabilities                        |
| `reports/summary.txt`                            | Output paths (also printed at completion)              |




## Train your own models

```bash
# 1) Compare assembly to the reference (once per assembly)
gffcompare -r reference.gtf -o cmp assembly.gtf
# produces cmp.assembly.gtf.tmap (name depends on gffcompare version / -o)

telos train \
  --bam sample.sorted.bam \
  --gtf assembly.gtf \
  --ref-gtf reference.gtf \
  --tmap cmp.assembly.gtf.tmap \
  --outdir train_out \
  --stage1-workers 8 \
  -v
```

Use the resulting `train_out/models/` directory with `telos predict`.

To train on **pooled assemblers** (recommended for paper-style models), concatenate GTFs/tmaps
with unique ID prefixes first, then pass the pooled files to `telos train`. 

### Training outputs


| Path                                          | Contents                                  |
| --------------------------------------------- | ----------------------------------------- |
| `models/…`                                    | Complete 8-file model bundle (see above)  |
| `predictions/sites.scored.tsv`                | Site scores on the training assembly      |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | Transcript rankings (+ diagnostic labels) |
| `reports/train_metrics.csv`                   | Stage I / Stage II validation metrics     |




## Pretrained models

We ship **code** in this repo and **model weights** separately (large joblibs; not in git).

```text
telos-models-<modality>/
  models/                 # the 8-file bundle
  config.yaml             # Stage I config used for training
  manifest.yaml           # train sample, assemblers, annotation version
```

Modalities trained for v2 (pooled assemblers, GENCODE v49):


| Modality      | Train sample | Assemblers           | Use for          |
| ------------- | ------------ | -------------------- | ---------------- |
| `ont-cdna`    | ENCFF023EXJ  | StringTie + IsoQuant | ONT cDNA         |
| `ont-drna`    | NA12878      | StringTie + IsoQuant | ONT direct RNA   |
| `pacbio`      | ENCFF450VAU  | StringTie + IsoQuant | PacBio Iso-Seq   |
| `short-reads` | SRR307903    | StringTie + Scallop2 | Illumina RNA-seq |


After download:

```bash
telos predict \
  --bam your.sorted.bam \
  --gtf your.assembly.gtf \
  --model-dir telos-models-ont-cdna/models \
  --outdir predict_out
```

Pick the modality that best matches your library type. Prefer training on matched data when
your protocol or species differs substantially.

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

