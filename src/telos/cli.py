"""
Command-line entry for Telos.

Defines ``argparse`` subcommands and maps parsed arguments to frozen IO dataclasses
(:mod:`telos.config_models`) and the thin handlers in :mod:`telos.commands`.
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

from telos import __version__
from telos.commands.benchmark import run_benchmark
from telos.commands.benchmark_matrix import run_benchmark_matrix
from telos.commands.predict import run_predict
from telos.commands.train import run_train
from telos.config_loader import default_stage1_config_path
from telos.config_models import BenchmarkIO, PredictIO, TrainIO
from telos.logging_setup import configure_logging


def _default_outdir(command: str) -> Path:
    """
    Build a timestamped output directory under the current working directory.

    Pattern: ``./telos_<command>_<YYYYMMDD_HHMMSS>/`` (resolved absolute). Used when ``train``,
    ``predict``, or ``benchmark`` omits ``--outdir``.
    """
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return (Path.cwd() / f"telos_{command}_{stamp}").resolve()


def build_parser() -> argparse.ArgumentParser:
    """
    Construct the root parser with four required subcommands.

    Subcommands:

    - ``train`` — paths for BAM, assembly GTF, reference GTF, gffcompare tmap; optional config and Stage I parallelism flags.
    - ``predict`` — BAM, GTF, ``--model-dir``; optional config and Stage I parallelism flags.
    - ``benchmark`` — benchmark YAML path; optional outdir.
    - ``benchmark-matrix`` — data type and train/test annotation choices; required benchmark ``--outdir``; optional bundle root and stage1 config.

    Returns:
        Configured :class:`argparse.ArgumentParser` (caller must invoke ``parse_args``).
    """
    p = argparse.ArgumentParser(
        prog="telos",
        description="Telos: train models and score assembled transcripts.",
    )
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(
        dest="command",
        required=True,
        metavar="{train,predict}",
    )

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument(
        "-v",
        "--verbose",
        action="count",
        default=0,
        help="Increase diagnostic detail; repeat for debug logging.",
    )
    common.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress non-error diagnostics (final output paths are still shown).",
    )

    train = sub.add_parser("train", parents=[common], help="Train Stage I/II models")
    train.add_argument("--bam", type=Path, required=True)
    train.add_argument("--gtf", type=Path, required=True)
    train.add_argument("--ref-gtf", type=Path, required=True)
    train.add_argument(
        "--tmap",
        type=Path,
        required=True,
        help="gffcompare .tmap (qry_id / class_code) for Stage II transcript labels.",
    )
    train.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: ./telos_train_<timestamp>).",
    )
    train.add_argument(
        "--config",
        type=Path,
        default=None,
        help="YAML/JSON config (default: bundled configs/stage1.defaults.yaml).",
    )
    train.add_argument(
        "--stage1-no-parallel",
        action="store_true",
        help="Disable multiprocessing for Stage I BAM feature extraction.",
    )
    train.add_argument(
        "--stage1-workers",
        type=int,
        default=None,
        metavar="N",
        help="Stage I feature pool size (default: config or min(CPU, 8)). Implies parallel when >1.",
    )
    predict = sub.add_parser(
        "predict", parents=[common], help="Run inference with trained models"
    )
    predict.add_argument("--bam", type=Path, required=True)
    predict.add_argument("--gtf", type=Path, required=True)
    predict.add_argument("--model-dir", type=Path, required=True)
    predict.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: ./telos_predict_<timestamp>).",
    )
    predict.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Optional YAML/JSON config for feature/model parameters.",
    )
    predict.add_argument(
        "--stage1-no-parallel",
        action="store_true",
        help="Disable multiprocessing for Stage I BAM feature extraction.",
    )
    predict.add_argument(
        "--stage1-workers",
        type=int,
        default=None,
        metavar="N",
        help="Stage I feature pool size (default: config or min(CPU, 8)). Implies parallel when >1.",
    )
    predict.add_argument(
        "--backend",
        choices=("xgb", "rf", "both"),
        default="xgb",
        help="Backend(s) for user-facing scored/filtered GTF output (default: xgb).",
    )
    predict.add_argument(
        "--min-score",
        type=float,
        default=None,
        help="Also write a GTF containing transcripts at or above this score.",
    )
    bench = sub.add_parser("benchmark", help=argparse.SUPPRESS)
    bench.add_argument("--config", type=Path, required=True)
    bench.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help="Output directory (default: ./telos_benchmark_<timestamp>).",
    )

    bmat = sub.add_parser(
        "benchmark-matrix",
        help=argparse.SUPPRESS,
    )
    bmat.add_argument(
        "--data-type",
        choices=("sr", "cdna", "drna", "pacbio"),
        required=True,
        help="Sequencing modality (fixed train sample per type).",
    )
    bmat.add_argument(
        "--train-annotation",
        choices=("refseq", "gencode", "ensembl"),
        required=True,
        dest="train_annotation",
        help="Reference bundle for training (BAM/GTF/tmap/ref GTF).",
    )
    bmat.add_argument(
        "--test-annotation",
        choices=("refseq", "gencode", "ensembl"),
        required=True,
        dest="test_annotation",
        help="Reference bundle for test samples (same modality, other samples).",
    )
    bmat.add_argument("--outdir", type=Path, required=True, help="Benchmark root (train/, tests/, reports/).")
    bmat.add_argument(
        "--bundles-root",
        type=Path,
        default=None,
        help="Directory containing <ref_id>/<modality>/<sample>/ (default: TELOS_BUNDLES_ROOT or ./data/bundles).",
    )
    bmat.add_argument(
        "--stage1-config",
        type=Path,
        default=None,
        help="Telos Stage I YAML (default: bundled telos/configs/stage1.defaults.yaml).",
    )

    # argparse.SUPPRESS does not remove subparser pseudo-actions on all supported
    # Python versions, so prune only their help entries while keeping parsing intact.
    sub._choices_actions = [
        action
        for action in sub._choices_actions
        if action.dest not in {"benchmark", "benchmark-matrix"}
    ]

    return p


def main(argv: list[str] | None = None) -> int:
    """
    Parse CLI arguments and dispatch to the appropriate command handler.

    Exit codes are delegated to ``run_train``, ``run_predict``, ``run_benchmark``, or
    ``run_benchmark_matrix`` (typically ``0`` success, ``1`` partial benchmark failure, ``2`` config/preflight).

    Args:
        argv: Argument vector; ``None`` means use ``sys.argv[1:]`` inside :func:`argparse.ArgumentParser.parse_args`.

    Returns:
        Integer exit code suitable for :func:`sys.exit`.
    """
    args = build_parser().parse_args(argv)
    configure_logging(args)
    if args.command == "train":
        outdir = args.outdir if args.outdir is not None else _default_outdir("train")
        train_cfg = args.config if args.config is not None else default_stage1_config_path()
        return run_train(
            TrainIO(
                bam=args.bam,
                gtf=args.gtf,
                ref_gtf=args.ref_gtf,
                tmap=args.tmap,
                outdir=outdir,
                config_file=train_cfg,
                stage1_no_parallel=args.stage1_no_parallel,
                stage1_n_workers=args.stage1_workers,
            )
        )
    if args.command == "predict":
        outdir = args.outdir if args.outdir is not None else _default_outdir("predict")
        return run_predict(
            PredictIO(
                bam=args.bam,
                gtf=args.gtf,
                model_dir=args.model_dir,
                outdir=outdir,
                config_file=args.config,
                backend=args.backend,
                min_score=args.min_score,
                stage1_no_parallel=args.stage1_no_parallel,
                stage1_n_workers=args.stage1_workers,
            )
        )
    if args.command == "benchmark":
        outdir = args.outdir if args.outdir is not None else _default_outdir("benchmark")
        return run_benchmark(BenchmarkIO(config=args.config, outdir=outdir))
    if args.command == "benchmark-matrix":
        return run_benchmark_matrix(
            data_type=args.data_type,
            train_annotation=args.train_annotation,
            test_annotation=args.test_annotation,
            outdir=args.outdir,
            bundles_root=args.bundles_root,
            stage1_config=args.stage1_config,
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
