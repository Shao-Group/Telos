# Telos

Telos scores transcript start sites (TSS), transcript end sites (TES), and complete
transcripts from an aligned RNA-seq BAM and an assembled GTF. It supports assemblies from
any tool that writes transcript-level `cov` attributes.

## Install

```bash
conda install -c bioconda telos
```

To install the current checkout:

```bash
python -m pip install .
```

Telos requires Python 3.10 or newer. Training and prediction are implemented in Python and
do not require RNASeqtools or gffcompare.

## Inputs

- A coordinate-sorted BAM with a `.bai` or `.csi` index.
- An assembled GTF whose transcript records contain numeric `cov` attributes.
- For training, a reference annotation GTF and the assembly's gffcompare `.tmap`.
- For prediction, a model directory produced by `telos train`. It must contain both RF and
  XGBoost Stage I bundles and both Stage II bundles.

## Train

```bash
telos train \
  --bam sample.bam \
  --gtf assembly.gtf \
  --ref-gtf reference.gtf \
  --tmap assembly.gtf.tmap \
  --outdir train_run
```

Useful options:

- `--config FILE`: override the bundled Stage I defaults.
- `--stage1-workers N` or `--stage1-no-parallel`: control BAM feature extraction.
- `--save-intermediates`: write optional diagnostics under `debug/`.
- `-v` / `-vv`: show progress or debug detail; `--quiet` suppresses non-error diagnostics.

Training always fits RF and XGBoost Stage I models and one Stage II model driven by each
backend.

### Training outputs

| Path | Contents |
| --- | --- |
| `models/stage1_{tss,tes}_{rf,xgb}_model.joblib` | Stage I model bundles |
| `models/stage2_model_{rf,xgb}.joblib` | Stage II models |
| `models/stage2_feature_names_{rf,xgb}.json` | Stage II feature order |
| `predictions/sites.scored.tsv` | Site identity columns plus `p_site_rf`, `p_site_xgb` |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | `transcript_id`, prediction, and diagnostic true label |
| `reports/train_metrics.csv` | Stage I and Stage II validation metrics |
| `reports/run_manifest.json` | Arguments, versions, input fingerprints, and split settings |

## Predict

```bash
telos predict \
  --bam sample.bam \
  --gtf assembly.gtf \
  --model-dir train_run/models \
  --backend xgb \
  --outdir predict_run
```

Use the same `--config` used during training whenever feature-extraction settings were
customized.

`--backend {xgb,rf,both}` selects the backend used for user-facing GTF output and defaults
to `xgb`. Prediction still scores both backends and writes both ranked TSVs. Add
`--min-score FLOAT` to write a high-confidence filtered GTF.

### Prediction outputs

| Path | Contents |
| --- | --- |
| `predictions/transcripts.scored.<backend>.gtf` | Input GTF with transcript `cov` replaced by Telos score |
| `predictions/transcripts.filtered.<backend>.gtf` | Score-filtered GTF, when `--min-score` is set |
| `predictions/transcripts.ranked.{rf,xgb}.tsv` | Both backend transcript rankings |
| `predictions/sites.scored.tsv` | Site scores with both backend probability columns |
| `reports/run_manifest.json` | Reproducibility metadata |
| `reports/summary.txt` | Concise output-path summary, also printed at completion |

With `--backend both`, Telos writes `.rf.gtf` and `.xgb.gtf` variants.

## Benchmarking and reproduction

Benchmarking, plotting, and paper reproduction workflows live in the
[Telos-test repository](https://github.com/Shao-Group/Telos-test). The legacy
`benchmark` and `benchmark-matrix` command names remain as redirects for compatibility.

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
