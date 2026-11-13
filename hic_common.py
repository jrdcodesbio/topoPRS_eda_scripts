"""Readers for the ENCODE/Juicer Hi-C call files (loops and contact-domain BEDPE).

BEDPE layout (observed in the QC run): the first line is the header, starting with
'#chr1'; later lines starting with '#' are '# juicer_tools ...' comments. Coordinates
are zero-based, half-open. Domain files repeat the name 'score', so duplicate names
are made unique ('score', 'score_2') rather than inferred from a data row.
"""
from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

from gwas_common import normalize_chrom

BEDPE_COORDS = ["chr1", "x1", "x2", "chr2", "y1", "y2"]


def unique_names(names: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out = []
    for n in names:
        seen[n] = seen.get(n, 0) + 1
        out.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    return out


def read_bedpe(path: Path) -> tuple[pd.DataFrame, dict]:
    """Return (table, metadata). Adds int64 coordinates and normalized chromosomes '_c1'/'_c2'."""
    opener = gzip.open if str(path).endswith(".gz") else open
    comments, rows = [], []
    with opener(path, "rt") as fh:
        header = fh.readline().rstrip("\n")
        if not header.startswith("#chr1"):
            raise ValueError(f"{path}: first line should be the '#chr1 ...' header, got {header[:80]!r}")
        names = unique_names(header.lstrip("#").split("\t"))
        if names[:6] != BEDPE_COORDS:
            raise ValueError(f"{path}: first six columns are {names[:6]}, expected {BEDPE_COORDS}")
        for line in fh:
            if line.startswith("#"):
                comments.append(line.rstrip("\n"))
            elif line.strip():
                rows.append(line.rstrip("\n").split("\t"))
    bad = [i for i, r in enumerate(rows) if len(r) != len(names)]
    if bad:
        raise ValueError(f"{path}: {len(bad)} rows have a column count different from the header ({len(names)})")
    df = pd.DataFrame(rows, columns=names)
    for c in ("x1", "x2", "y1", "y2"):
        df[c] = pd.to_numeric(df[c], errors="raise").astype("int64")
    df["_c1"] = normalize_chrom(df["chr1"])
    df["_c2"] = normalize_chrom(df["chr2"])
    meta = {"path": str(path), "n_rows": len(df), "columns": names, "comments": comments}
    return df, meta
