#!/usr/bin/env python3
"""Match significant GWAS variants to the LD reference .pvar and write the plink2 --clump input.

A GWAS variant matches the reference when chromosome (with 'chr' stripped), position, and
the unordered allele pair are all the same. The clump file uses the *reference* ID, so
plink2 can find each variant. All original GWAS columns are kept.

Outputs:
  <out-clump>                       ID, P, REF_CHROM, REF_POS + all original GWAS columns
  <qc-dir>/unmatched_sig_snps.tsv   significant GWAS rows not found in the reference
                                    (plink2 drops these, so they're listed here instead)
  <qc-dir>/dropped_incomplete.tsv   significant rows missing chrom/pos/allele fields
  <qc-dir>/prep_summary.json
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import pandas as pd

from gwas_common import (
    DEFAULT_COLS,
    allele_pair_key,
    is_missing,
    normalize_chrom,
    read_gwas_chunks,
    read_pvar_chunks,
    require_columns,
    setup_logging,
)

log = logging.getLogger("prep_clump_input")

OUT_KEY_COLS = ["ID", "P", "REF_CHROM", "REF_POS"]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gwas", type=Path, required=True)
    p.add_argument("--ref-pvar", type=Path, required=True, help="Uncompressed .pvar of the LD reference")
    p.add_argument("--out-clump", type=Path, required=True)
    p.add_argument("--qc-dir", type=Path, required=True)
    p.add_argument("--p-max", type=float, required=True,
                   help="Keep GWAS rows with P <= this (should be max of --clump-p1/--clump-p2)")
    p.add_argument("--chrom-col", default=DEFAULT_COLS["chrom"])
    p.add_argument("--pos-col", default=DEFAULT_COLS["pos"])
    p.add_argument("--ea-col", default=DEFAULT_COLS["effect_allele"])
    p.add_argument("--oa-col", default=DEFAULT_COLS["other_allele"])
    p.add_argument("--p-col", default=DEFAULT_COLS["p"])
    p.add_argument("--chunksize", type=int, default=1_000_000)
    return p.parse_args()


def select_significant(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Return (significant rows with match keys, significant rows that are incomplete, total rows read)."""
    keep, incomplete, n_rows = [], [], 0
    for chunk in read_gwas_chunks(args.gwas, args.chunksize):
        if n_rows == 0:
            require_columns(chunk, [args.chrom_col, args.pos_col, args.ea_col, args.oa_col, args.p_col],
                            str(args.gwas))
            clash = [c for c in OUT_KEY_COLS if c in chunk.columns]
            if clash:
                raise ValueError(f"GWAS columns {clash} collide with output key columns {OUT_KEY_COLS}")
        n_rows += len(chunk)
        p = pd.to_numeric(chunk[args.p_col].where(~is_missing(chunk[args.p_col])), errors="coerce")
        sig = chunk[p <= args.p_max].copy()
        if sig.empty:
            continue
        sig["P"] = p[p <= args.p_max]
        bad = (is_missing(sig[args.chrom_col]) | is_missing(sig[args.ea_col]) | is_missing(sig[args.oa_col])
               | pd.to_numeric(sig[args.pos_col], errors="coerce").isna())
        incomplete.append(sig[bad])
        sig = sig[~bad]
        sig["_chrom"] = normalize_chrom(sig[args.chrom_col])
        sig["_pos"] = pd.to_numeric(sig[args.pos_col]).astype("int64")
        sig["_alleles"] = allele_pair_key(sig[args.ea_col], sig[args.oa_col])
        keep.append(sig)
    if n_rows == 0:
        raise ValueError(f"{args.gwas}: file is empty")
    sig = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    bad = pd.concat(incomplete, ignore_index=True) if incomplete else pd.DataFrame()
    return sig, bad, n_rows


def lookup_reference(pvar: Path, wanted: pd.DataFrame, chunksize: int) -> pd.DataFrame:
    """Return the reference variants at the wanted (chrom, pos) positions, with match keys."""
    wanted_pos = set(zip(wanted["_chrom"], wanted["_pos"]))
    hits = []
    for chunk in read_pvar_chunks(pvar, chunksize):
        chrom = normalize_chrom(chunk["CHROM"])
        pos = pd.to_numeric(chunk["POS"]).astype("int64")
        mask = pd.Series([k in wanted_pos for k in zip(chrom, pos)], index=chunk.index)
        if mask.any():
            h = chunk[mask].copy()
            h["_chrom"] = chrom[mask]
            h["_pos"] = pos[mask]
            # Multiallelic records ("A,T") can't be matched as a pair; they are left out
            # here and therefore appear in the unmatched list.
            h = h[~h["ALT"].str.contains(",")]
            h["_alleles"] = allele_pair_key(h["REF"], h["ALT"])
            hits.append(h)
    if not hits:
        return pd.DataFrame(columns=["CHROM", "POS", "ID", "REF", "ALT", "_chrom", "_pos", "_alleles"])
    return pd.concat(hits, ignore_index=True)


def check_clump_ids(pvar: Path, ids: pd.Series, chunksize: int, qc_dir: Path) -> None:
    """plink2 --clump looks variants up by ID, so matched IDs must be real and unique in the reference."""
    n_blank = int(ids.isin([".", ""]).sum())
    if n_blank:
        raise SystemExit(f"{n_blank} matched reference variants have ID '.' or empty; plink2 --clump "
                         "cannot use them. Reset IDs in the reference first.")
    wanted = set(ids)
    counts = pd.Series(dtype="int64")
    for chunk in read_pvar_chunks(pvar, chunksize):
        counts = counts.add(chunk.loc[chunk["ID"].isin(wanted), "ID"].value_counts(), fill_value=0)
    dup = counts[counts > 1]
    if len(dup):
        path = qc_dir / "duplicate_ref_ids.tsv"
        dup.astype("int64").rename("n_records").to_csv(path, sep="\t", header=True)
        raise SystemExit(f"{len(dup)} matched reference IDs occur more than once in {pvar}; see {path}. "
                         "plink2 --clump needs unique IDs.")


def match(sig: pd.DataFrame, ref: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, int]:
    """Join on (chrom, pos, alleles). Return (clump input, unmatched rows, number of duplicate rows dropped)."""
    keys = ["_chrom", "_pos", "_alleles"]
    ref = ref.drop_duplicates(keys + ["ID"])
    merged = sig.merge(ref[keys + ["ID", "CHROM", "POS"]], on=keys, how="left", indicator=True)
    unmatched = merged[merged["_merge"] == "left_only"].drop(columns=["ID", "CHROM", "POS", "_merge"])
    matched = merged[merged["_merge"] == "both"].drop(columns=["_merge"])
    # plink2 needs one row per ID. Several GWAS rows on the same reference variant keep the smallest P.
    n_before = len(matched)
    matched = matched.sort_values("P", kind="stable").drop_duplicates("ID", keep="first")
    n_dup = n_before - len(matched)
    matched = matched.rename(columns={"CHROM": "REF_CHROM", "POS": "REF_POS"})
    gwas_cols = [c for c in sig.columns if c not in keys and c != "P"]
    out = matched[OUT_KEY_COLS + gwas_cols].sort_values(["REF_CHROM", "REF_POS"], kind="stable")
    return out, unmatched.drop(columns=keys), n_dup


def main() -> None:
    setup_logging()
    args = parse_args()
    args.qc_dir.mkdir(parents=True, exist_ok=True)
    args.out_clump.parent.mkdir(parents=True, exist_ok=True)

    sig, incomplete, n_rows = select_significant(args)
    log.info("read %d GWAS rows; %d at P <= %g (%d incomplete)", n_rows, len(sig) + len(incomplete),
             args.p_max, len(incomplete))
    if sig.empty:
        raise SystemExit(f"No complete GWAS rows with P <= {args.p_max:g}; nothing to clump.")

    ref = lookup_reference(args.ref_pvar, sig, args.chunksize)
    out, unmatched, n_dup = match(sig, ref)
    check_clump_ids(args.ref_pvar, out["ID"], args.chunksize, args.qc_dir)
    frac_matched = (len(sig) - len(unmatched)) / len(sig)
    if frac_matched < 0.5:
        log.warning("only %.1f%% of significant variants matched the reference; check the genome build "
                    "and variant coverage of the reference", 100 * frac_matched)

    # 10 significant digits keeps near-tied P values in the same order for plink2.
    out.to_csv(args.out_clump, sep="\t", index=False, float_format="%.10g")
    unmatched.to_csv(args.qc_dir / "unmatched_sig_snps.tsv", sep="\t", index=False, float_format="%.10g")
    incomplete.to_csv(args.qc_dir / "dropped_incomplete.tsv", sep="\t", index=False, float_format="%.10g")

    summary = {
        "gwas": str(args.gwas),
        "ref_pvar": str(args.ref_pvar),
        "p_max": args.p_max,
        "n_gwas_rows": n_rows,
        "n_sig_complete": len(sig),
        "n_sig_incomplete": len(incomplete),
        "n_ref_records_at_sig_positions": len(ref),
        "n_matched_unique_ids": len(out),
        "n_duplicate_rows_collapsed": n_dup,
        "n_unmatched": len(unmatched),
        "frac_sig_matched": frac_matched,
        "unmatched_by_chrom": unmatched[args.chrom_col].value_counts().to_dict() if len(unmatched) else {},
    }
    with open(args.qc_dir / "prep_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    log.info("matched %d unique IDs, %d unmatched (%.1f%% matched) -> %s", len(out), len(unmatched),
             100 * frac_matched, args.out_clump)


if __name__ == "__main__":
    main()
