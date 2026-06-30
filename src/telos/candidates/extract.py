"""
Extract TSS and TES **candidate sites** from transcript lines in an assembly GTF.

Each transcript yields two :class:`CandidateSite` rows (one TSS, one TES) using the same 1-based
coordinate convention as rnaseqtools/gtfformat TSSTES on forward and reverse strands.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from telos.gtf_attributes import parse_transcript_id


@dataclass(frozen=True)
class CandidateSite:
    """
    One genomic site (TSS or TES) belonging to a transcript.

    ``position`` is 1-based on the reference, consistent with GTF ``start``/``end`` parsing below.
    ``strand`` is ``+`` or ``-`` as in column 7 of the source GTF.
    """

    transcript_id: str
    site_type: str  # TSS or TES
    chrom: str
    position: int
    strand: str


def extract_candidate_sites_from_gtf(gtf_path: Path) -> list[CandidateSite]:
    """
    Scan ``gtf_path`` for ``feature == transcript`` rows and emit two candidates per transcript.

    Coordinate rules (1-based):

    - ``+`` strand: TSS at ``start + 1``, TES at ``end`` (inclusive end as in GTF).
    - ``-`` strand: TSS at ``end``, TES at ``start + 1``.

    Strand handling:

    - ``+`` strand: emit TSS at ``start + 1`` and TES at ``end``.
    - ``-`` strand: emit TSS at ``end`` and TES at ``start + 1``.
    - ``.`` strand (unknown): emit **both** sets above (two strands × TSS/TES), so downstream
      joins never silently drop transcripts with missing strand annotations.

    Skips lines with fewer than 9 columns, non-transcript features, or missing ``transcript_id``.

    Returns:
        List in file order (TSS then TES for each kept transcript).
    """
    candidates: list[CandidateSite] = []
    with gtf_path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attrs = cols
            if feature != "transcript":
                continue
            tx_id = parse_transcript_id(attrs)
            if not tx_id:
                continue

            s = int(start)
            e = int(end)

            def _emit_for_canonical_strand(canon_strand: str) -> None:
                if canon_strand == "+":
                    tss_pos = s + 1
                    tes_pos = e
                else:
                    tss_pos = e
                    tes_pos = s + 1

                candidates.append(
                    CandidateSite(
                        transcript_id=tx_id,
                        site_type="TSS",
                        chrom=chrom,
                        position=tss_pos,
                        strand=canon_strand,
                    )
                )
                candidates.append(
                    CandidateSite(
                        transcript_id=tx_id,
                        site_type="TES",
                        chrom=chrom,
                        position=tes_pos,
                        strand=canon_strand,
                    )
                )

            if strand in {"+", "-"}:
                _emit_for_canonical_strand(strand)
            elif strand == ".":
                # Ambiguous strand: keep the transcript for both coordinate conventions.
                _emit_for_canonical_strand("+")
                _emit_for_canonical_strand("-")
            else:
                # Unknown symbol (rare): skip.
                continue
    return candidates


def write_candidates_tsv(candidates: list[CandidateSite], out_tsv: Path) -> None:
    """
    Write a tab-separated table of candidates with header ``transcript_id, site_type, chrom, position, strand``.

    Creates parent directories. Used for debugging or external tool handoff.
    """
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["transcript_id", "site_type", "chrom", "position", "strand"])
        for c in candidates:
            writer.writerow([c.transcript_id, c.site_type, c.chrom, c.position, c.strand])
