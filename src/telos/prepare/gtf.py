"""
Normalize assembler GTFs for Telos and optionally build a gffcompare ``.tmap``.

Telos expects each transcript feature to carry a numeric ``cov "..."`` attribute.
IsoQuant often stores abundance in a separate TPM TSV (and sometimes uses
``coverage`` instead of ``cov``). Training also needs a gffcompare ``.tmap`` for
Stage II labels; this module shells out to ``gffcompare`` when requested.
"""

from __future__ import annotations

import logging
import re
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from telos.gtf_attributes import parse_transcript_id

logger = logging.getLogger(__name__)

_COV_RE = re.compile(r'(?<![A-Za-z0-9_])cov "([^"]*)"')
_COVERAGE_RE = re.compile(r'(?<![A-Za-z0-9_])coverage "([^"]*)"')


class PrepareGtfError(RuntimeError):
    """Raised when GTF preparation or gffcompare fails."""


@dataclass(frozen=True)
class PrepareGtfResult:
    """Paths and counts produced by :func:`prepare_assembly_gtf`."""

    out_gtf: Path
    out_tmap: Path | None
    transcript_lines: int
    updated_cov: int
    missing_abundance: int


def load_abundance_table(path: Path) -> dict[str, float]:
    """
    Parse a two-column transcript abundance TSV (IsoQuant TPM style).

    Accepts optional ``#`` comment / header lines. Column 1 = transcript id,
    column 2 = numeric abundance written into ``cov``.
    """
    path = path.resolve()
    if not path.is_file():
        raise PrepareGtfError(f"Abundance table not found: {path}")

    mapping: dict[str, float] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            txid = parts[0].strip()
            if not txid:
                continue
            try:
                value = float(parts[1].strip())
            except ValueError:
                continue
            mapping[txid] = value
    if not mapping:
        raise PrepareGtfError(f"No transcript abundance rows parsed from: {path}")
    return mapping


def _set_cov_attr(attributes: str, value: float) -> str:
    """Write ``cov "<value>"``, replacing ``cov`` or ``coverage`` if present."""
    cov_literal = f'cov "{value:.6f}"'
    if _COV_RE.search(attributes):
        return _COV_RE.sub(cov_literal, attributes, count=1)
    if _COVERAGE_RE.search(attributes):
        return _COVERAGE_RE.sub(cov_literal, attributes, count=1)
    attrs = attributes.rstrip()
    if attrs and not attrs.endswith(";"):
        attrs = f"{attrs};"
    if attrs:
        return f"{attrs} {cov_literal};"
    return f"{cov_literal};"


def _normalize_coverage_key(attributes: str) -> tuple[str, bool]:
    """
    Rename ``coverage "..."`` → ``cov "..."`` when ``cov`` is absent.

    Returns ``(attributes, changed)``.
    """
    if _COV_RE.search(attributes):
        return attributes, False
    match = _COVERAGE_RE.search(attributes)
    if not match:
        return attributes, False
    value = match.group(1)
    try:
        float(value)
    except ValueError as exc:
        raise PrepareGtfError(f"Non-numeric coverage attribute: {attributes!r}") from exc
    # Replace key name while keeping the quoted value.
    new_attrs = _COVERAGE_RE.sub(f'cov "{value}"', attributes, count=1)
    return new_attrs, True


def inject_cov_from_tpm(
    gtf_path: Path,
    tpm_path: Path,
    out_gtf_path: Path,
) -> PrepareGtfResult:
    """Rewrite transcript ``cov`` from an IsoQuant-style TPM / abundance TSV."""
    abundance = load_abundance_table(tpm_path)
    return _rewrite_transcript_cov(
        gtf_path=gtf_path,
        out_gtf_path=out_gtf_path,
        abundance=abundance,
        normalize_coverage_key=True,
    )


def normalize_transcript_cov_attrs(
    gtf_path: Path,
    out_gtf_path: Path,
) -> PrepareGtfResult:
    """
    Copy ``gtf_path`` → ``out_gtf_path``, renaming ``coverage`` → ``cov`` when needed.

    Does not invent abundance values; transcripts without ``cov``/``coverage`` are
    left unchanged (Telos preflight / cov parsing will fail later if still missing).
    """
    return _rewrite_transcript_cov(
        gtf_path=gtf_path,
        out_gtf_path=out_gtf_path,
        abundance=None,
        normalize_coverage_key=True,
    )


def _rewrite_transcript_cov(
    *,
    gtf_path: Path,
    out_gtf_path: Path,
    abundance: dict[str, float] | None,
    normalize_coverage_key: bool,
) -> PrepareGtfResult:
    gtf_path = gtf_path.resolve()
    out_gtf_path = out_gtf_path.resolve()
    if not gtf_path.is_file():
        raise PrepareGtfError(f"GTF not found: {gtf_path}")

    out_gtf_path.parent.mkdir(parents=True, exist_ok=True)
    same_path = out_gtf_path == gtf_path

    transcript_lines = 0
    updated_cov = 0
    missing_abundance = 0

    final_out = out_gtf_path
    temp_out: Path | None = None
    if same_path:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=str(out_gtf_path.parent),
            delete=False,
            prefix=f"{out_gtf_path.name}.tmp.",
        ) as handle:
            temp_out = Path(handle.name)
        final_out = temp_out

    with gtf_path.open("r", encoding="utf-8", errors="replace") as src, final_out.open(
        "w", encoding="utf-8"
    ) as dst:
        for raw in src:
            line = raw.rstrip("\n")
            if not line or line.startswith("#"):
                dst.write(raw if raw.endswith("\n") else raw + "\n")
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "transcript":
                dst.write(raw if raw.endswith("\n") else raw + "\n")
                continue

            transcript_lines += 1
            attrs = cols[8]
            changed = False

            if abundance is not None:
                txid = parse_transcript_id(attrs)
                if txid is None:
                    dst.write("\t".join(cols) + "\n")
                    continue
                value = abundance.get(txid)
                if value is None:
                    missing_abundance += 1
                    if normalize_coverage_key:
                        attrs, key_changed = _normalize_coverage_key(attrs)
                        changed = key_changed
                else:
                    attrs = _set_cov_attr(attrs, value)
                    changed = True
                    updated_cov += 1
            elif normalize_coverage_key:
                attrs, key_changed = _normalize_coverage_key(attrs)
                if key_changed:
                    changed = True
                    updated_cov += 1

            if changed:
                cols[8] = attrs
                dst.write("\t".join(cols) + "\n")
            else:
                dst.write(raw if raw.endswith("\n") else raw + "\n")

    if same_path and temp_out is not None:
        temp_out.replace(out_gtf_path)

    return PrepareGtfResult(
        out_gtf=out_gtf_path,
        out_tmap=None,
        transcript_lines=transcript_lines,
        updated_cov=updated_cov,
        missing_abundance=missing_abundance,
    )


def _resolve_gffcompare(gffcompare_bin: Path | None) -> Path:
    if gffcompare_bin is not None:
        path = gffcompare_bin.resolve()
        if not path.is_file():
            raise PrepareGtfError(f"gffcompare not found: {path}")
        return path
    found = shutil.which("gffcompare")
    if not found:
        raise PrepareGtfError(
            "gffcompare not found on PATH. Install it (e.g. "
            "`conda install -c bioconda gffcompare`) or pass --gffcompare /path/to/gffcompare."
        )
    return Path(found)


def run_gffcompare_tmap(
    *,
    query_gtf: Path,
    ref_gtf: Path,
    outdir: Path,
    gffcompare_bin: Path | None = None,
    tmap_name: str = "assembly.tmap",
) -> Path:
    """
    Run gffcompare and copy the resulting ``.tmap`` to ``outdir / tmap_name``.

    Uses local basenames inside ``outdir`` so gffcompare 0.11.x writes a ``.tmap``
    (absolute query paths can skip the tmap file).
    """
    query_gtf = query_gtf.resolve()
    ref_gtf = ref_gtf.resolve()
    outdir = outdir.resolve()
    if not query_gtf.is_file():
        raise PrepareGtfError(f"Query GTF not found: {query_gtf}")
    if not ref_gtf.is_file():
        raise PrepareGtfError(f"Reference GTF not found: {ref_gtf}")

    gffcompare = _resolve_gffcompare(gffcompare_bin)
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / "gffcompare_work"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    query_link = work / "query.gtf"
    ref_link = work / "ref.gtf"
    query_link.symlink_to(query_gtf)
    ref_link.symlink_to(ref_gtf)

    prefix = "cmp"
    cmd = [
        str(gffcompare),
        "-r",
        "ref.gtf",
        "-o",
        prefix,
        "query.gtf",
    ]
    logger.info("Running: %s (cwd=%s)", " ".join(cmd), work)
    try:
        completed = subprocess.run(
            cmd,
            cwd=str(work),
            check=False,
            capture_output=True,
            text=True,
        )
    except OSError as exc:
        raise PrepareGtfError(f"Failed to execute gffcompare: {exc}") from exc

    if completed.returncode != 0:
        detail = (completed.stderr or completed.stdout or "").strip()
        raise PrepareGtfError(
            f"gffcompare failed (exit {completed.returncode})"
            + (f": {detail}" if detail else "")
        )

    tmaps = sorted(work.glob("*.tmap"))
    if not tmaps:
        raise PrepareGtfError(
            f"gffcompare finished but no .tmap was written under {work}. "
            "Check that gffcompare supports -r/-o with relative query paths."
        )
    # Prefer the query.tmap style when multiple sidecars exist.
    preferred = [p for p in tmaps if "query" in p.name] or tmaps
    src_tmap = preferred[0]
    dest = outdir / tmap_name
    shutil.copy2(src_tmap, dest)
    logger.info("Wrote tmap: %s (from %s)", dest, src_tmap.name)
    return dest


def prepare_assembly_gtf(
    *,
    gtf: Path,
    outdir: Path,
    tpm: Path | None = None,
    ref_gtf: Path | None = None,
    make_tmap: bool = False,
    gffcompare_bin: Path | None = None,
    out_gtf_name: str = "assembly.gtf",
    out_tmap_name: str = "assembly.tmap",
) -> PrepareGtfResult:
    """
    Prepare an assembly GTF for Telos and optionally build ``assembly.tmap``.

    - With ``tpm``: inject abundance into transcript ``cov``.
    - Without ``tpm``: copy GTF and rename ``coverage`` → ``cov`` when present.
    - With ``make_tmap``: require ``ref_gtf`` and run gffcompare.
    """
    outdir = outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    out_gtf = outdir / out_gtf_name

    if tpm is not None:
        result = inject_cov_from_tpm(gtf, tpm, out_gtf)
    else:
        result = normalize_transcript_cov_attrs(gtf, out_gtf)

    out_tmap: Path | None = None
    if make_tmap:
        if ref_gtf is None:
            raise PrepareGtfError("--make-tmap requires --ref-gtf")
        out_tmap = run_gffcompare_tmap(
            query_gtf=out_gtf,
            ref_gtf=ref_gtf,
            outdir=outdir,
            gffcompare_bin=gffcompare_bin,
            tmap_name=out_tmap_name,
        )

    return PrepareGtfResult(
        out_gtf=result.out_gtf,
        out_tmap=out_tmap,
        transcript_lines=result.transcript_lines,
        updated_cov=result.updated_cov,
        missing_abundance=result.missing_abundance,
    )
