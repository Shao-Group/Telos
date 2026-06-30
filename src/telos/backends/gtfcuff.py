"""
Python implementation of gtfcuff-style ROC/AUC sweeps over gffcompare ``.tmap`` rows.

"""

from __future__ import annotations

from pathlib import Path


class GtfcuffError(RuntimeError):
    """Raised when a tmap file is missing or unparsable for gtfcuff-style metrics."""


def _read_tmap_items(tmap_path: Path) -> list[dict[str, object]]:
    """Parse a gffcompare ``.tmap`` into dict rows with class code, abundance fields, and length."""
    tmap_path = tmap_path.resolve()
    if not tmap_path.is_file():
        raise GtfcuffError(f"TMAP not found: {tmap_path}")
    out: list[dict[str, object]] = []
    with tmap_path.open("r", encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 12:
                parts = line.split()
            if len(parts) < 12:
                continue
            if parts[0] == "ref_gene_id":
                continue
            try:
                out.append(
                    {
                        "ref_transcript_id": parts[1],
                        "code": (parts[2] or "@")[0],
                        "fpkm": float(parts[6]),
                        "tpm": float(parts[7]),
                        "coverage": float(parts[8]),
                        "length": int(float(parts[9])),
                    }
                )
            except ValueError:
                continue
    return out


def run_gtfcuff_roc(
    gtfcuff_bin: Path | None,
    tmap_path: Path,
    ref_multi_exon_count: int,
    measure: str = "cov",
) -> str:
    """Sweep sorted tmap rows and emit recall/precision TSV text (header + numeric rows)."""
    _ = gtfcuff_bin
    refsize = int(ref_multi_exon_count)
    items = _read_tmap_items(tmap_path)
    if not items:
        return ""

    m = (measure or "cov").strip()
    if m == "TPM":
        key = "tpm"
    elif m == "FPKM":
        key = "fpkm"
    else:
        key = "coverage"
    vt = sorted(items, key=lambda it: float(it[key]))

    mt: dict[str, int] = {}
    for it in vt:
        if str(it["code"]) == "=":
            rt = str(it["ref_transcript_id"])
            mt[rt] = mt.get(rt, 0) + 1

    correct = len(mt)
    total_correct = correct
    change = True
    lines: list[str] = ["recall\tprecision"]
    n = len(vt)
    for i, it in enumerate(vt):
        rem = n - i
        sen_pct = (correct * 100.0 / refsize) if refsize > 0 else 0.0
        pre_pct = (correct * 100.0 / rem) if rem > 0 else 0.0
        if ((total_correct - correct) % 10 == 0) and change:
            # Emit a point on the same step cadence as the original gtfcuff tool for comparability.
            lines.append(f"{sen_pct / 100.0:.6f}\t{pre_pct / 100.0:.6f}")

        if str(it["code"]) == "=":
            rt = str(it["ref_transcript_id"])
            c = mt.get(rt, 0)
            if c <= 1:
                if rt in mt:
                    del mt[rt]
                correct -= 1
            else:
                mt[rt] = c - 1
            change = True
        else:
            change = False
    return "\n".join(lines) + "\n"


def run_gtfcuff_auc(gtfcuff_bin: Path | None, tmap_path: Path, ref_multi_exon_count: int) -> float:
    """Compute the unnormalized AUC statistic used by the original gtfcuff ``auc`` mode."""
    _ = gtfcuff_bin
    refsize = int(ref_multi_exon_count)
    items = sorted(_read_tmap_items(tmap_path), key=lambda it: float(it["coverage"]))
    if not items:
        return 0.0

    correct = sum(1 for it in items if str(it["code"]) == "=")
    sen0 = (correct * 100.0 / refsize) if refsize > 0 else 0.0
    pre0 = (correct * 100.0 / len(items)) if items else 0.0
    auc = sen0 * pre0

    n = len(items)
    for i in range(n - 1):
        if str(items[i]["code"]) == "=":
            correct -= 1
        den = n - i - 1
        sen = (correct * 100.0 / refsize) if refsize > 0 else 0.0
        pre = (correct * 100.0 / den) if den > 0 else 0.0
        area = (sen + sen0) * 0.5 * (pre - pre0)
        auc += area
        pre0 = pre
        sen0 = sen
    return float(auc)
