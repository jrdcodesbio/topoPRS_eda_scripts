#!/usr/bin/env python3
"""Schema and distribution QC for the GWAS summary statistics (chunked; makes no changes).

Writes to --out-dir:
  column_summary.tsv   per-column missing counts and the raw missing tokens seen
  value_counts.tsv     value counts for low-cardinality fields (chromosome, hm_* codes)
  summary.json         row counts, P-value range, allele checks, duplicates, significance counts
  head.tsv             first rows, as read
"""
from __future__ import annotations

import argparse
import json
import logging
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from gwas_common import (
    DEFAULT_COLS,
    MISSING_TOKENS,
    allele_pair_key,
    is_missing,
    normalize_chrom,
    read_gwas_chunks,
    require_columns,
    setup_logging,
)

log = logging.getLogger("inspect_gwas")

VALUE_COUNT_COLS = ["hm_code", "hm_coordinate_conversion"]  # the --chrom-col column is added at runtime
ACGT = r"^[ACGT]$"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gwas", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--chrom-col", default=DEFAULT_COLS["chrom"])
    p.add_argument("--pos-col", default=DEFAULT_COLS["pos"])
    p.add_argument("--ea-col", default=DEFAULT_COLS["effect_allele"])
    p.add_argument("--oa-col", default=DEFAULT_COLS["other_allele"])
    p.add_argument("--p-col", default=DEFAULT_COLS["p"])
    p.add_argument("--thresholds", type=float, nargs="+", default=[5e-8, 1e-6, 1e-5])
    p.add_argument("--chunksize", type=int, default=1_000_000)
    p.add_argument("--head-rows", type=int, default=20)
    return p.parse_args()


def inspect(args: argparse.Namespace) -> dict:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    n_rows = 0
    missing = Counter()
    missing_tokens: dict[str, Counter] = {}
    value_counts: dict[str, Counter] = {}
    p_stats = {"n_numeric": 0, "n_unparseable": 0, "n_zero": 0, "n_gt1": 0, "n_neg": 0,
               "min_positive": np.inf, "max": -np.inf}
    sig_counts = Counter()
    allele_stats = Counter()
    pos_keys: list[pd.DataFrame] = []
    columns: list[str] | None = None

    for i, chunk in enumerate(read_gwas_chunks(args.gwas, args.chunksize)):
        if columns is None:
            columns = list(chunk.columns)
            require_columns(chunk, [args.chrom_col, args.pos_col, args.ea_col, args.oa_col, args.p_col],
                            str(args.gwas))
            chunk.head(args.head_rows).to_csv(args.out_dir / "head.tsv", sep="\t", index=False)
            missing_tokens = {c: Counter() for c in columns}
            value_counts = {c: Counter() for c in dict.fromkeys([args.chrom_col, *VALUE_COUNT_COLS]) if c in columns}
        n_rows += len(chunk)

        for c in columns:
            stripped = chunk[c].str.strip()
            m = stripped.isin(MISSING_TOKENS)
            missing[c] += int(m.sum())
            missing_tokens[c].update(stripped[m].tolist())
        for c in value_counts:
            value_counts[c].update(chunk[c].tolist())

        # P values: parse numerically, but count anything unparseable rather than drop it silently.
        p_raw = chunk[args.p_col]
        p_present = ~is_missing(p_raw)
        p = pd.to_numeric(p_raw.where(p_present), errors="coerce")
        p_stats["n_unparseable"] += int((p_present & p.isna()).sum())
        pv = p.dropna()
        p_stats["n_numeric"] += len(pv)
        p_stats["n_zero"] += int((pv == 0).sum())
        p_stats["n_gt1"] += int((pv > 1).sum())
        p_stats["n_neg"] += int((pv < 0).sum())
        if (pv > 0).any():
            p_stats["min_positive"] = min(p_stats["min_positive"], float(pv[pv > 0].min()))
        if len(pv):
            p_stats["max"] = max(p_stats["max"], float(pv.max()))
        for t in args.thresholds:
            sig_counts[f"p<{t:g}"] += int((pv < t).sum())

        ea = chunk[args.ea_col].str.strip().str.upper()
        oa = chunk[args.oa_col].str.strip().str.upper()
        snv = ea.str.match(ACGT) & oa.str.match(ACGT)
        allele_stats["biallelic_snv"] += int(snv.sum())
        allele_stats["non_snv_or_nonacgt"] += int((~snv).sum())
        allele_stats["ea_equals_oa"] += int((ea == oa).sum())
        allele_stats["palindromic_AT_CG"] += int((snv & (allele_pair_key(ea, oa).isin(["A:T", "C:G"]))).sum())

        keys = pd.DataFrame({
            "chrom": normalize_chrom(chunk[args.chrom_col]),
            "pos": pd.to_numeric(chunk[args.pos_col], errors="coerce"),
            "alleles": allele_pair_key(ea, oa),
        })
        allele_stats["pos_unparseable"] += int(keys["pos"].isna().sum())
        pos_keys.append(keys)
        log.info("chunk %d: %d rows (total %d)", i, len(chunk), n_rows)

    if columns is None:
        raise ValueError(f"{args.gwas}: file is empty")

    keys = pd.concat(pos_keys, ignore_index=True)
    dup = {
        "dup_chrom_pos_rows": int(keys.duplicated(["chrom", "pos"], keep=False).sum()),
        "dup_chrom_pos_alleles_rows": int(keys.duplicated(["chrom", "pos", "alleles"], keep=False).sum()),
    }
    chrom_norm_counts = keys["chrom"].value_counts().to_dict()

    pd.DataFrame({
        "column": columns,
        "n_missing": [missing[c] for c in columns],
        "frac_missing": [missing[c] / n_rows for c in columns],
        "missing_tokens": [json.dumps(dict(missing_tokens[c])) for c in columns],
    }).to_csv(args.out_dir / "column_summary.tsv", sep="\t", index=False)

    vc_rows = [{"column": c, "value": v, "count": n}
               for c, cnt in value_counts.items() for v, n in sorted(cnt.items(), key=lambda x: -x[1])]
    pd.DataFrame(vc_rows, columns=["column", "value", "count"]).to_csv(
        args.out_dir / "value_counts.tsv", sep="\t", index=False)

    for k in ("min_positive", "max"):
        if not np.isfinite(p_stats[k]):
            p_stats[k] = None
    summary = {
        "gwas": str(args.gwas),
        "n_rows": n_rows,
        "columns": columns,
        "p_value": p_stats,
        "significance_counts": dict(sig_counts),
        "alleles": dict(allele_stats),
        "duplicates": dup,
        "normalized_chrom_counts": chrom_norm_counts,
    }
    with open(args.out_dir / "summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    return summary


def main() -> None:
    setup_logging()
    args = parse_args()
    summary = inspect(args)
    log.info("rows=%d sig=%s p_zero=%d -> %s", summary["n_rows"], summary["significance_counts"],
             summary["p_value"]["n_zero"], args.out_dir)


if __name__ == "__main__":
    main()
