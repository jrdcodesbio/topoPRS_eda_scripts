"""Shared helpers for reading the GWAS summary statistics and plink2 .pvar files.

Every field is read as a string and missing tokens are handled explicitly. This
keeps the QC able to report how missing values are actually encoded, and stops
rsIDs, alleles, or chromosome codes from being coerced to other types.
"""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Iterator

import pandas as pd

# Tokens counted as missing. QC also reports how often each one appears raw.
MISSING_TOKENS = frozenset({"", "NA", "NaN", "nan", "#NA", ".", "None"})

# Defaults match the header of sourceData/GWAS/allen3wayGWAS.h.tsv.gz.
DEFAULT_COLS = {
    "chrom": "chromosome",
    "pos": "base_pair_location",
    "effect_allele": "effect_allele",
    "other_allele": "other_allele",
    "p": "p_value",
}

log = logging.getLogger(__name__)


def setup_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )


def read_gwas_chunks(path: Path, chunksize: int = 1_000_000) -> Iterator[pd.DataFrame]:
    """Yield the tab-separated GWAS file as all-string chunks with no NA conversion."""
    return pd.read_csv(
        path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        chunksize=chunksize,
        compression="infer",
    )


def is_missing(s: pd.Series) -> pd.Series:
    return s.str.strip().isin(MISSING_TOKENS)


def normalize_chrom(s: pd.Series) -> pd.Series:
    """Map 'chr1'/'1'/'Chr1' to '1' and 'chrX'/'23' to 'X', so the GWAS and the reference can be matched."""
    out = s.str.strip().str.replace(r"^(?i:chr)", "", regex=True).str.upper()
    return out.replace({"23": "X", "24": "Y", "25": "XY", "26": "MT", "M": "MT"})


def allele_pair_key(a: pd.Series, b: pd.Series) -> pd.Series:
    """Build an order-independent allele key, e.g. ('G','A') -> 'A:G'.

    LD clumping only needs to identify the variant, not the effect direction, so
    REF/ALT vs effect/other ordering doesn't matter. Strand flips are not
    resolved: non-palindromic SNPs on the opposite strand will fail to match and
    are reported as unmatched.
    """
    a = a.str.strip().str.upper()
    b = b.str.strip().str.upper()
    lo = a.where(a <= b, b)
    hi = b.where(a <= b, a)
    return lo + ":" + hi


def read_pvar_header(path: Path) -> tuple[list[str], int]:
    """Return the .pvar column names and the number of lines up to and including the '#CHROM' header."""
    with open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                return line.lstrip("#").rstrip("\n").split("\t"), i + 1
            raise ValueError(f"{path}: expected '#CHROM' header line, found: {line[:80]!r}")
    raise ValueError(f"{path}: no '#CHROM' header line found")


def read_pvar_chunks(path: Path, chunksize: int = 2_000_000) -> Iterator[pd.DataFrame]:
    """Yield CHROM/POS/ID/REF/ALT chunks from an uncompressed plink2 .pvar."""
    cols, skip = read_pvar_header(path)
    needed = ["CHROM", "POS", "ID", "REF", "ALT"]
    missing = [c for c in needed if c not in cols]
    if missing:
        raise ValueError(f"{path}: missing pvar columns {missing}; header={cols}")
    return pd.read_csv(
        path,
        sep="\t",
        header=None,
        names=cols,
        usecols=needed,
        skiprows=skip,
        dtype=str,
        keep_default_na=False,
        na_filter=False,
        chunksize=chunksize,
    )


def require_columns(df: pd.DataFrame, cols: list[str], source: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise KeyError(f"{source}: missing required columns {missing}; available={list(df.columns)}")
