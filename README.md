# Telos

Telos scores transcript start sites (TSS), transcript end sites (TES), and assembled
transcripts from a coordinate-sorted BAM and a transcript GTF. It works with any
assembler that provides transcript abundance (StringTie, IsoQuant, Scallop2, …).

```text
prepare-gtf  →  (optional) normalize cov / build tmap
train        →  models/
predict      →  scored & ranked transcripts
```

## Install

```bash
git clone https://github.com/Shao-Group/Telos.git
cd Telos
conda env create -f environment.yml
conda activate telos
pip install -e . --no-deps
telos --version
```

Requires **Python ≥ 3.10**. The conda env also installs `gffcompare` (used by
`telos prepare-gtf --make-tmap`).

When published on Bioconda:

```bash
conda install -c bioconda -c conda-forge telos
```



## Commands at a glance


| Command             | Purpose                                                          |
| ------------------- | ---------------------------------------------------------------- |
| `telos prepare-gtf` | Make an assembly Telos-ready (`cov` attribute; optional `.tmap`) |
| `telos train`       | Train Stage I + Stage II models                                  |
| `telos predict`     | Score a new BAM/GTF with a trained model directory               |




### Inputs


| Command         | Required                                                                                |
| --------------- | --------------------------------------------------------------------------------------- |
| **predict**     | Sorted BAM + index (`.bai`/`.csi`), assembly GTF with transcript `cov "..."`, `models/` |
| **train**       | Same BAM/GTF, plus reference GTF and gffcompare `.tmap`                                 |
| **prepare-gtf** | Assembler GTF; optional IsoQuant TPM; optional `--make-tmap` + reference GTF            |


A complete model directory has 8 files:

```text
models/
  stage1_tss_{rf,xgb}_model.joblib
  stage1_tes_{rf,xgb}_model.joblib
  stage2_model_{rf,xgb}.joblib
  stage2_feature_names_{rf,xgb}.json
```



## 1. Prepare the assembly GTF

Telos reads abundance from transcript attribute `cov "..."`.


| Assembler            | What to do                                               |
| -------------------- | -------------------------------------------------------- |
| StringTie / Scallop2 | Usually already have `cov`. For training, build a tmap.  |
| IsoQuant             | Inject TPM into `cov`, then (for training) build a tmap. |


```bash
# IsoQuant: TPM → cov
telos prepare-gtf \
  --gtf OUT.transcript_models.gtf \
  --tpm OUT.transcript_model_tpm.tsv \
  --outdir prepared

# StringTie / Scallop2: build tmap for training
telos prepare-gtf \
  --gtf stringtie.gtf \
  --ref-gtf gencode.gtf \
  --make-tmap \
  --outdir prepared

# IsoQuant: TPM → cov and tmap
telos prepare-gtf \
  --gtf OUT.transcript_models.gtf \
  --tpm OUT.transcript_model_tpm.tsv \
  --ref-gtf gencode.gtf \
  --make-tmap \
  --outdir prepared
# → prepared/assembly.gtf  and  prepared/assembly.tmap
```



## 2. Predict (pretrained or your own models)

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

Extra Options:

- `--backend {xgb,rf,both}` — which backend writes the scored GTF (default `xgb`)
- `--min-score 0.8` — also write a filtered high-confidence GTF containing transcripts crossing the min-score.
- `--config FILE` — use the **same** Stage I feature config as training if you customized it
- `--stage1-workers N` — BAM feature-extraction parallelism



### Prediction outputs


| Path                                             | Contents                                |
| ------------------------------------------------ | --------------------------------------- |
| `predictions/transcripts.scored.<backend>.gtf`   | GTF with transcript `cov` = Telos score |
| `predictions/transcripts.filtered.<backend>.gtf` | Only if `--min-score` is set            |
| `predictions/transcripts.ranked.{rf,xgb}.tsv`    | Ranked scores                           |
| `predictions/sites.scored.tsv`                   | Site-level probabilities                |
| `reports/summary.txt`                            | Output path summary                     |




## 3. Train your own models

```bash
telos train \
  --bam sample.sorted.bam \
  --gtf prepared/assembly.gtf \
  --ref-gtf gencode.gtf \
  --tmap prepared/assembly.tmap \
  --outdir train_out \
  --split-policy chr1-10 \
  --n-jobs -1 \
  --stage1-workers 8 \
  -v
```

Then run predict with `--model-dir train_out/models`.

### Important training options


| Flag               | Default          | Meaning                                                                                                                                                |
| ------------------ | ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `--split-policy`   | `chr1-10`        | Primary autosomes in this range = **train**; all other contigs (remaining autosomes, X/Y/MT, scaffolds) = **validation**. Accepts `chr1-10` or `1-10`. |
| `--n-jobs`         | config (`-1`)    | Threads for RF, XGBoost, and LightGBM (`-1` = all CPUs)                                                                                                |
| `--stage1-workers` | config           | Parallelism for BAM feature extraction                                                                                                                 |
| `--config`         | bundled defaults | Full Stage I YAML (hyperparameters, feature windows, …)                                                                                                |


Chromosome names are normalized across UCSC (`chr1`), Ensembl (`1`), and RefSeq (`NC_000001.11`).

### Advanced config

Copy and edit the packaged defaults, then pass `--config`:

```bash
# defaults live at: $(python -c "import telos, pathlib; print(pathlib.Path(telos.__file__).parent/'configs'/'stage1.defaults.yaml')")
telos train ... --config my_stage1.yaml
```

Use the **same** `--config` for `predict` if you changed feature-extraction settings. Model hyperparameters and window sizes are YAML-only.

### Training outputs


| Path                                          | Contents                                 |
| --------------------------------------------- | ---------------------------------------- |
| `models/…`                                    | 8-file model bundle                      |
| `predictions/sites.scored.tsv`                | Site scores on the training assembly     |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | Ranked transcripts (+ diagnostic labels) |
| `reports/train_metrics.csv`                   | Stage I / Stage II validation metrics    |




## Pretrained models

Weights are **not** in git (large joblibs). Download a modality bundle from the GitHub
Release / Zenodo (links TBD) and point `--model-dir` at its `models/` folder.


| Bundle        | Train sample | Assemblers           | Use for          |
| ------------- | ------------ | -------------------- | ---------------- |
| `ont-cdna`    | ENCFF023EXJ  | StringTie + IsoQuant | ONT cDNA         |
| `ont-drna`    | NA12878      | StringTie + IsoQuant | ONT direct RNA   |
| `pacbio`      | ENCFF450VAU  | StringTie + IsoQuant | PacBio Iso-Seq   |
| `short-reads` | SRR307903    | StringTie + Scallop2 | Illumina RNA-seq |


```bash
telos predict \
  --bam your.sorted.bam \
  --gtf your.assembly.gtf \
  --model-dir telos-models-ont-cdna/models \
  --outdir predict_out
```

Prefer matched modality; retrain when species or protocol differ substantially.

Paper benchmarks: [Telos-test](https://github.com/Shao-Group/Telos-test).

## License and citation

BSD 3-Clause.

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

