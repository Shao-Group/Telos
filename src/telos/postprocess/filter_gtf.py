"""
Filter GTF features to high-confidence transcripts using Stage II ``transcripts.ranked*.tsv`` scores.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

_TX_RE = re.compile(r'transcript_id\s+"([^"]+)"')
_GENE_RE = re.compile(r'gene_id\s+"([^"]+)"')


def _extract_attr(attrs: str, pattern: re.Pattern[str]) -> str | None:
    """Return first regex capture from GTF column 9 ``attrs`` or ``None``."""
    m = pattern.search(attrs)
    return m.group(1) if m else None


def _load_kept_transcripts(ranked_tsv: Path, prob_threshold: float) -> set[str]:
    """
    ``transcript_id`` set where ``pred_prob >= prob_threshold`` in ranked TSV.

    Raises:
        ValueError: Missing required columns.
    """
    df = pd.read_csv(ranked_tsv, sep="\t")
    if "transcript_id" not in df.columns or "pred_prob" not in df.columns:
        raise ValueError("ranked TSV must contain transcript_id and pred_prob columns")
    kept = df[df["pred_prob"].astype(float) >= float(prob_threshold)]["transcript_id"].astype(str)
    return set(kept.tolist())


def filter_gtf_by_transcript_scores(
    in_gtf: Path,
    ranked_tsv: Path,
    out_gtf: Path,
    prob_threshold: float,
) -> tuple[int, int]:
    """
    Two-pass GTF filter: keep transcript features in ``keep_tx`` and any exon/CDS sharing their ids.

    First pass counts transcripts and collects ``gene_id`` for kept transcripts. Second pass writes
    lines: all genes that still have a kept transcript, all features with ``transcript_id`` in ``keep_tx``,
    and any feature tied to a kept ``gene_id`` (e.g. exons without transcript_id in attrs — rare).

    Args:
        in_gtf: Input GTF path.
        ranked_tsv: Stage II ranked output with ``transcript_id`` and ``pred_prob``.
        out_gtf: Output path (parent dirs created).
        prob_threshold: Minimum ``pred_prob`` to retain a transcript.

    Returns:
        ``(kept_transcript_count, total_transcript_count_seen)`` from the first pass.
    """
    keep_tx = _load_kept_transcripts(ranked_tsv, prob_threshold)
    if not in_gtf.is_file():
        raise FileNotFoundError(f"Input GTF not found: {in_gtf}")

    total_tx = 0
    kept_tx = 0
    keep_gene_ids: set[str] = set()

    with in_gtf.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw or raw.startswith("#"):
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            feat = cols[2]
            attrs = cols[8]
            tx = _extract_attr(attrs, _TX_RE)
            if feat == "transcript" and tx is not None:
                total_tx += 1
                if tx in keep_tx:
                    kept_tx += 1
                    gid = _extract_attr(attrs, _GENE_RE)
                    if gid:
                        keep_gene_ids.add(gid)

    out_gtf.parent.mkdir(parents=True, exist_ok=True)
    with in_gtf.open("r", encoding="utf-8", errors="replace") as src, out_gtf.open(
        "w", encoding="utf-8"
    ) as dst:
        for raw in src:
            if not raw:
                continue
            if raw.startswith("#"):
                dst.write(raw)
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9:
                dst.write(raw)
                continue
            feat = cols[2]
            attrs = cols[8]
            tx = _extract_attr(attrs, _TX_RE)
            gid = _extract_attr(attrs, _GENE_RE)

            if feat == "gene":
                if gid is None or gid in keep_gene_ids:
                    dst.write(raw)
                continue
            if tx is not None:
                if tx in keep_tx:
                    dst.write(raw)
                continue
            if gid is not None and gid in keep_gene_ids:
                dst.write(raw)

    return kept_tx, total_tx
