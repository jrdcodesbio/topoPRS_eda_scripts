#!/usr/bin/env python3
"""Convert the plink2 .clumps report into tidy locus tables for the Hi-C annotation step.

Outputs (in --out-dir):
  loci.tsv            one row per clump: lead variant, P, member count, member bp span
  locus_members.tsv   one row per (locus, variant), with the original GWAS columns joined in
Outputs (in --qc-dir):
  unmatched_near_loci.tsv  significant GWAS SNPs missing from the LD reference, with the
                           nearest locus lead and whether they fall within --kb of it
  clump_summary.json
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd

from gwas_common import DEFAULT_COLS, normalize_chrom, setup_logging

log = logging.getLogger("summarize_clumps")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clumps", type=Path, required=True, help="plink2 <prefix>.clumps")
    p.add_argument("--clump-input", type=Path, required=True, help="Output of prep_clump_input.py")
    p.add_argument("--unmatched", type=Path, required=True, help="unmatched_sig_snps.tsv from prep_clump_input.py")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--qc-dir", type=Path, required=True)
    p.add_argument("--kb", type=float, required=True, help="Clump window used, for flagging nearby unmatched SNPs")
    p.add_argument("--chrom-col", default=DEFAULT_COLS["chrom"])
    p.add_argument("--pos-col", default=DEFAULT_COLS["pos"])
    return p.parse_args()


def read_clumps(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    df.columns = [c.lstrip("#") for c in df.columns]
    missing = [c for c in ("CHROM", "POS", "ID", "P", "SP2") if c not in df.columns]
    if missing:
        raise ValueError(f"{path}: missing expected .clumps columns {missing}; header={list(df.columns)}")
    return df


def parse_sp2(field: str) -> list[str]:
    """Split an SP2 entry into variant IDs. '.' means no members; a '(...)' suffix carries file/A1 info and is dropped."""
    if field in ("", "."):
        return []
    return [tok.split("(", 1)[0] for tok in field.split(",") if tok]


def build_members(clumps: pd.DataFrame, clump_input: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for i, r in enumerate(clumps.itertuples(index=False), start=1):
        locus_id = f"L{i:04d}"
        rows.append((locus_id, r.ID, True))
        rows.extend((locus_id, m, False) for m in parse_sp2(r.SP2))
    members = pd.DataFrame(rows, columns=["locus_id", "ID", "is_lead"])
    out = members.merge(clump_input, on="ID", how="left", validate="many_to_one")
    n_missing = int(out["REF_POS"].isna().sum())
    if n_missing:
        raise ValueError(f"{n_missing} clump member IDs not found in clump input; inputs are out of sync")
    return out


def build_loci(members: pd.DataFrame) -> pd.DataFrame:
    m = members.assign(REF_POS=members["REF_POS"].astype("int64"), P=members["P"].astype(float))
    lead = m[m["is_lead"]].set_index("locus_id")
    agg = m.groupby("locus_id").agg(
        n_variants=("ID", "size"),
        span_start=("REF_POS", "min"),
        span_end=("REF_POS", "max"),
        member_ids=("ID", lambda s: ",".join(s)),
    )
    loci = lead[["ID", "REF_CHROM", "REF_POS", "P"]].rename(
        columns={"ID": "lead_id", "REF_CHROM": "chrom", "REF_POS": "lead_pos", "P": "lead_p"}).join(agg)
    if "rsid" in lead.columns:
        loci.insert(1, "lead_rsid", lead["rsid"])
    loci["span_bp"] = loci["span_end"] - loci["span_start"] + 1
    return loci.reset_index()


def flag_unmatched(unmatched: pd.DataFrame, loci: pd.DataFrame, kb: float,
                   chrom_col: str, pos_col: str) -> pd.DataFrame:
    """Annotate each reference-absent significant SNP with its nearest lead.

    These SNPs couldn't be clumped. One within the window of an existing lead may
    belong to that locus. One far from every lead may be a signal the LD reference
    misses entirely.
    """
    if unmatched.empty:
        return unmatched.assign(nearest_locus=pd.Series(dtype=str), dist_to_lead_bp=pd.Series(dtype=float),
                                within_kb_of_lead=pd.Series(dtype=bool))
    u = unmatched.copy()
    u["_chrom"] = normalize_chrom(u[chrom_col])
    u["_pos"] = pd.to_numeric(u[pos_col]).astype("int64")
    l_chrom = normalize_chrom(loci["chrom"].astype(str))
    nearest, dist = [], []
    for c, pos in zip(u["_chrom"], u["_pos"]):
        sub = loci[l_chrom == c]
        if sub.empty:
            nearest.append(None)
            dist.append(np.nan)
            continue
        d = (sub["lead_pos"] - pos).abs()
        j = d.idxmin()
        nearest.append(sub.at[j, "locus_id"])
        dist.append(float(d[j]))
    u["nearest_locus"] = nearest
    u["dist_to_lead_bp"] = dist
    u["within_kb_of_lead"] = u["dist_to_lead_bp"] <= kb * 1000
    return u.drop(columns=["_chrom", "_pos"])


def main() -> None:
    setup_logging()
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.qc_dir.mkdir(parents=True, exist_ok=True)

    clumps = read_clumps(args.clumps)
    clump_input = pd.read_csv(args.clump_input, sep="\t", dtype=str, keep_default_na=False)
    unmatched = pd.read_csv(args.unmatched, sep="\t", dtype=str, keep_default_na=False)

    members = build_members(clumps, clump_input)
    loci = build_loci(members)
    flagged = flag_unmatched(unmatched, loci, args.kb, args.chrom_col, args.pos_col)

    loci.to_csv(args.out_dir / "loci.tsv", sep="\t", index=False)
    members.to_csv(args.out_dir / "locus_members.tsv", sep="\t", index=False)
    flagged.to_csv(args.qc_dir / "unmatched_near_loci.tsv", sep="\t", index=False)

    n_in = len(clump_input)
    summary = {
        "n_loci": len(loci),
        "n_clump_input_variants": n_in,
        "n_variants_assigned_to_loci": len(members),
        "n_input_variants_not_in_any_clump": n_in - members["ID"].nunique(),
        "variants_per_locus": loci["n_variants"].describe().to_dict(),
        "span_bp": loci["span_bp"].describe().to_dict(),
        "n_unmatched_sig": len(flagged),
        "n_unmatched_within_kb_of_lead": int(flagged["within_kb_of_lead"].sum()) if len(flagged) else 0,
        "n_unmatched_far_from_all_leads": int((~flagged["within_kb_of_lead"].astype(bool)).sum()) if len(flagged) else 0,
    }
    with open(args.qc_dir / "clump_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2, default=float)
    log.info("%d loci from %d variants; %d reference-absent sig SNPs (%d far from any lead)",
             len(loci), n_in, len(flagged), summary["n_unmatched_far_from_all_leads"])


if __name__ == "__main__":
    main()
