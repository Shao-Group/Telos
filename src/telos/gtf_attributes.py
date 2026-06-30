"""
GTF column-9 attribute parsing shared across backends and evaluation.

Keeps one regex for ``transcript_id`` so BAM/GTF joins and gffcompare-related code agree on IDs.
"""

from __future__ import annotations

import re

# GTF9 attributes often encode ``transcript_id "TX.1";`` — only the double-quoted form is matched.
TRANSCRIPT_ID_RE = re.compile(r'transcript_id\s+"([^"]+)"')


def parse_transcript_id(attributes: str) -> str | None:
    """
    Extract the first ``transcript_id`` token from a GTF attributes string.

    Returns ``None`` when the pattern does not match (unquoted IDs, gene-only rows, etc.).
    """
    m = TRANSCRIPT_ID_RE.search(attributes)
    return m.group(1) if m else None
