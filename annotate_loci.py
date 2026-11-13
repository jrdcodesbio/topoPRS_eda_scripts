#!/usr/bin/env python3
"""Annotate clumped GWAS loci against Hi-C loops, contact domains and compartment scores.

Each file answers a different question, per reference (e.g. lung, Tcell):
  loops.bedpe.gz        does a lead/member variant sit on a loop anchor, and where does that contact lead?
                        (observed/expected/FDR describe the loop *call*, not enrichment of IPF variants)
  contactDom.bedpe.gz   which self-interacting domain(s) contain the variant, and how far is the nearest
                        domain boundary? (domains are diagonal intervals, not extra loop edges)
  compartments.bigWig   compartment score at the variant, labelled A (> 0), B (< 0), 'unassigned' (exactly 0),
                        or left missing (no coverage/NaN). Positive = A follows the orientation of the ENCODE
                        compartment tracks as stated by the user; the ENCODE pipeline docs reviewed do not state it.

Coordinates: variant POS is 1-based, so it occupies the 0-based base p0 = POS-1; BEDPE intervals are
0-based half-open, so a variant overlaps [start, end) when start <= p0 < end.

Outputs (--out-dir): variant_annotations.tsv, locus_loops.tsv, locus_domains.tsv, loci_annotated.tsv
QC (--qc-dir):       annotate_summary.json
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import re
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from gwas_common import normalize_chrom, setup_logging
from hic_common import read_bedpe

log = logging.getLogger("annotate_loci")

# (chrom as in the variant table, p0) -> (score, bin_start, bin_end); NaN/None when no value.
Lookup = Callable[[str, int], tuple]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--loci", type=Path, required=True, help="loci.tsv from summarize_clumps.py")
    p.add_argument("--members", type=Path, required=True, help="locus_members.tsv from summarize_clumps.py")
    p.add_argument("--reference", action="append", required=True, metavar="NAME=DIR",
                   help="Hi-C reference directory with loops.bedpe.gz, contactDom.bedpe.gz, compartments.bigWig "
                        "(repeatable)")
    p.add_argument("--near-kb", type=float, default=10, help="'near' an anchor if within this distance (no overlap)")
    p.add_argument("--loop-stat-regex", default=r"^(observed|expected.*|fdr.*)$",
                   help="Loop columns to carry into the output (matched against the actual header)")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--qc-dir", type=Path, required=True)
    return p.parse_args()


def load_variants(members: pd.DataFrame) -> pd.DataFrame:
    v = pd.DataFrame({
        "locus_id": members["locus_id"],
        "ID": members["ID"],
        "rsid": members["rsid"] if "rsid" in members.columns else "",
        "is_lead": members["is_lead"].astype(str).str.lower() == "true",
        "P": pd.to_numeric(members["P"]),
        "chrom": normalize_chrom(members["REF_CHROM"].astype(str)),
        "pos": pd.to_numeric(members["REF_POS"]).astype("int64"),
    })
    v["p0"] = v["pos"] - 1
    return v.reset_index(drop=True)


def interval_distance(p0: int, start: np.ndarray, end: np.ndarray) -> np.ndarray:
    """Bases between p0 and [start, end): 0 inside, start-p0 before, p0-(end-1) after."""
    return np.maximum(np.maximum(start - p0, p0 - (end - 1)), 0)


def loop_hits(variants: pd.DataFrame, loops: pd.DataFrame, near_bp: int) -> tuple[pd.DataFrame, pd.Series]:
    """Return (hit rows, spanning-loop count per variant).

    Hit rows: one per (variant, loop) with relation 'anchor' (overlaps anchor 1 and/or 2) or 'near'
    (within near_bp of an anchor without overlapping). 'Spanning' (strictly between the two anchors of a
    cis loop) is only counted, so it is never mistaken for contact.
    """
    c1, c2 = loops["_c1"].to_numpy(), loops["_c2"].to_numpy()
    x1, x2, y1, y2 = (loops[c].to_numpy() for c in ("x1", "x2", "y1", "y2"))
    rows, spanning = [], []
    for i, v in variants.iterrows():
        on1, on2 = c1 == v.chrom, c2 == v.chrom
        d1 = np.where(on1, interval_distance(v.p0, x1, x2), np.iinfo(np.int64).max)
        d2 = np.where(on2, interval_distance(v.p0, y1, y2), np.iinfo(np.int64).max)
        a1, a2 = d1 == 0, d2 == 0
        near = ~a1 & ~a2 & ((d1 <= near_bp) | (d2 <= near_bp))
        spanning.append(int((on1 & on2 & ~a1 & ~a2 & (v.p0 >= x2) & (v.p0 < y1)).sum()))
        for j in np.flatnonzero(a1 | a2):
            anchor = "both" if a1[j] and a2[j] else ("1" if a1[j] else "2")
            rows.append((i, j, "anchor", anchor, 0))
        for j in np.flatnonzero(near):
            anchor = "1" if d1[j] <= d2[j] else "2"
            rows.append((i, j, "near", anchor, int(min(d1[j], d2[j]))))
    hits = pd.DataFrame(rows, columns=["v_idx", "loop_idx", "relation", "anchor", "dist_to_anchor_bp"])
    return hits, pd.Series(spanning, index=variants.index, dtype="int64")


def partner_of(loops: pd.DataFrame, loop_idx: int, anchor: str) -> tuple[str, int, int]:
    """The other end of the loop. For 'both' (anchors overlap each other) anchor 2 is reported."""
    r = loops.iloc[loop_idx]
    if anchor == "2":
        return r["chr1"], int(r["x1"]), int(r["x2"])
    return r["chr2"], int(r["y1"]), int(r["y2"])


def merge_intervals(intervals: list[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    """Merge overlapping/abutting half-open intervals per chromosome."""
    out: list[list] = []
    for c, s, e in sorted(intervals, key=lambda t: (t[0], t[1], t[2])):
        if out and out[-1][0] == c and s <= out[-1][2]:
            out[-1][2] = max(out[-1][2], e)
        else:
            out.append([c, s, e])
    return [tuple(x) for x in out]


def domain_annotations(variants: pd.DataFrame, domains: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Per-variant domain fields, and (variant, domain) containment rows. Nested domains are all reported."""
    dc, s, e = domains["_c1"].to_numpy(), domains["x1"].to_numpy(), domains["x2"].to_numpy()
    per_var, contain = [], []
    for i, v in variants.iterrows():
        on = dc == v.chrom
        inside = np.flatnonzero(on & (s <= v.p0) & (v.p0 < e))
        bounds = np.concatenate([s[on], e[on]])
        if len(bounds):
            k = int(np.argmin(np.abs(bounds - v.p0)))
            nearest, dist = int(bounds[k]), int(abs(bounds[k] - v.p0))
        else:
            nearest, dist = None, None
        spans = e[inside] - s[inside]
        per_var.append({
            "n_domains": len(inside),
            "domain_ids": ",".join(f"D{j}" for j in inside),
            "smallest_domain_span_bp": int(spans.min()) if len(inside) else None,
            "nearest_domain_boundary": nearest,
            "dist_to_domain_boundary_bp": dist,
        })
        contain.extend((i, j) for j in inside)
    per_var = pd.DataFrame(per_var, index=variants.index)
    for c in ("smallest_domain_span_bp", "nearest_domain_boundary", "dist_to_domain_boundary_bp"):
        per_var[c] = per_var[c].astype("Int64")  # nullable: no domains on a chromosome -> <NA>, not a float
    return per_var, pd.DataFrame(contain, columns=["v_idx", "domain_idx"])


def compartment_annotations(variants: pd.DataFrame, lookup: Lookup) -> pd.DataFrame:
    out = []
    for _, v in variants.iterrows():
        score, b0, b1 = lookup(v.chrom, int(v.p0))
        has = score is not None and not math.isnan(score)
        out.append({
            "compartment_score": score if has else None,
            "compartment_bin_start": b0 if has else None,
            "compartment_bin_end": b1 if has else None,
            "compartment_sign": ("+" if score > 0 else "-" if score < 0 else "0") if has else "NA",
            "compartment_label": ("A" if score > 0 else "B" if score < 0 else "unassigned") if has else None,
        })
    out = pd.DataFrame(out, index=variants.index)
    for c in ("compartment_bin_start", "compartment_bin_end"):
        out[c] = out[c].astype("Int64")
    return out


def bigwig_lookup(path: Path) -> tuple[Lookup, dict]:
    import pyBigWig  # imported here so the rest of the module (and its tests) don't need it

    bw = pyBigWig.open(str(path))
    if bw is None or not bw.isBigWig():
        raise ValueError(f"{path}: not a readable bigWig")
    chroms = bw.chroms()
    name_map = {normalize_chrom(pd.Series([c])).iloc[0]: c for c in chroms}
    widths: set[int] = set()

    def lookup(chrom: str, p0: int):
        name = name_map.get(chrom)
        if name is None or p0 >= chroms[name]:
            return float("nan"), None, None
        iv = bw.intervals(name, p0, p0 + 1)
        if not iv:
            return float("nan"), None, None
        s, e, val = iv[0]
        widths.add(int(e - s))
        return float(val), int(s), int(e)

    meta = {"path": str(path), "n_chroms": len(chroms), "chrom_names_sample": list(chroms)[:5],
            "observed_bin_widths": widths}
    return lookup, meta


def annotate_reference(name: str, variants: pd.DataFrame, loops: pd.DataFrame, domains: pd.DataFrame,
                       lookup: Lookup, near_bp: int, stat_cols: list[str]):
    """Return (variant table, locus-loop table, locus-domain table) for one reference."""
    hits, spanning = loop_hits(variants, loops, near_bp)
    dom, contain = domain_annotations(variants, domains)
    comp = compartment_annotations(variants, lookup)

    anchor_hits = hits[hits["relation"] == "anchor"]
    var = variants.drop(columns="p0").copy()
    var.insert(0, "reference", name)
    var["n_anchor_loops"] = anchor_hits.groupby("v_idx").size().reindex(variants.index, fill_value=0)
    var["n_near_loops"] = hits[hits["relation"] == "near"].groupby("v_idx").size().reindex(
        variants.index, fill_value=0)
    var["n_spanning_loops"] = spanning
    var = pd.concat([var, dom, comp], axis=1)

    lead_idx = set(variants.index[variants["is_lead"]])
    loop_rows = []
    for (locus, j), g in hits.assign(locus_id=variants.loc[hits["v_idx"], "locus_id"].to_numpy()).groupby(
            ["locus_id", "loop_idx"]):
        on = g[g["relation"] == "anchor"]
        use = on if len(on) else g
        anchor = use["anchor"].iloc[0]
        pc, ps, pe = partner_of(loops, j, anchor)
        loop_rows.append({
            "reference": name, "locus_id": locus, "loop_id": f"loop{j}",
            "relation": "anchor" if len(on) else "near",
            "anchors_hit": ",".join(sorted(set(use["anchor"]))),
            # Both ends inside the locus: the contact stays within the locus, so the "partner" is not a
            # distal target.
            "both_anchors_in_locus": bool(len(on)) and ({"1", "2"} <= set(on["anchor"]) or "both" in set(on["anchor"])),
            "variant_ids": ",".join(variants.loc[use["v_idx"], "ID"]),
            "lead_on_anchor": bool(set(on["v_idx"]) & lead_idx),
            "min_dist_to_anchor_bp": int(use["dist_to_anchor_bp"].min()),
            "anchor1": f"{loops.iloc[j]['chr1']}:{loops.iloc[j]['x1']}-{loops.iloc[j]['x2']}",
            "anchor2": f"{loops.iloc[j]['chr2']}:{loops.iloc[j]['y1']}-{loops.iloc[j]['y2']}",
            "partner_chrom": pc, "partner_start": ps, "partner_end": pe,
            **{c: loops.iloc[j][c] for c in stat_cols},
        })
    loop_tab = pd.DataFrame(loop_rows)

    dom_rows = []
    if len(contain):
        c = contain.assign(locus_id=variants.loc[contain["v_idx"], "locus_id"].to_numpy())
        for (locus, j), g in c.groupby(["locus_id", "domain_idx"]):
            d = domains.iloc[j]
            dom_rows.append({
                "reference": name, "locus_id": locus, "domain_id": f"D{j}", "chrom": d["chr1"],
                "start": int(d["x1"]), "end": int(d["x2"]), "span_bp": int(d["x2"] - d["x1"]),
                "variant_ids": ",".join(variants.loc[g["v_idx"], "ID"]),
                "contains_lead": bool(variants.loc[g["v_idx"], "is_lead"].any()),
            })
    return var, loop_tab, pd.DataFrame(dom_rows)


def summarize_loci(loci: pd.DataFrame, var: pd.DataFrame, loop_tab: pd.DataFrame, name: str,
                   fdr_cols: list[str]) -> pd.DataFrame:
    rows = []
    for _, L in loci.iterrows():
        mv = var[var["locus_id"] == L["locus_id"]]
        lead = mv[mv["is_lead"]].iloc[0]
        lt = loop_tab[loop_tab["locus_id"] == L["locus_id"]] if len(loop_tab) else loop_tab
        on = lt[lt["relation"] == "anchor"] if len(lt) else lt
        partners = merge_intervals(list(zip(on["partner_chrom"], on["partner_start"], on["partner_end"]))) \
            if len(on) else []
        scores = pd.to_numeric(mv["compartment_score"], errors="coerce")
        row = {
            "reference": name, **L.to_dict(),
            "lead_on_anchor": bool(lead["n_anchor_loops"] > 0),
            "n_members_on_anchor": int((mv["n_anchor_loops"] > 0).sum()),
            "n_anchor_loops": len(on),
            "n_anchor_loops_both_ends_in_locus": int(on["both_anchors_in_locus"].sum()) if len(on) else 0,
            "n_near_only_loops": int((lt["relation"] == "near").sum()) if len(lt) else 0,
            "n_partner_regions": len(partners),
            "partner_regions": ";".join(f"{c}:{s}-{e}" for c, s, e in partners),
            "lead_n_domains": int(lead["n_domains"]),
            "lead_smallest_domain_span_bp": lead["smallest_domain_span_bp"],
            "lead_dist_to_domain_boundary_bp": lead["dist_to_domain_boundary_bp"],
            "lead_compartment_score": lead["compartment_score"],
            "lead_compartment_sign": lead["compartment_sign"],
            "members_compartment_min": scores.min(),
            "members_compartment_max": scores.max(),
            "lead_compartment_label": lead["compartment_label"],
            "n_members_A": int((mv["compartment_label"] == "A").sum()),
            "n_members_B": int((mv["compartment_label"] == "B").sum()),
            "n_members_unassigned": int((mv["compartment_label"] == "unassigned").sum()),
            "n_members_compartment_na": int(scores.isna().sum()),
        }
        for c in fdr_cols:
            row[f"min_{c}_anchor_loops"] = pd.to_numeric(on[c], errors="coerce").min() if len(on) else None
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    setup_logging()
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.qc_dir.mkdir(parents=True, exist_ok=True)
    near_bp = int(args.near_kb * 1000)

    members = pd.read_csv(args.members, sep="\t", dtype=str, keep_default_na=False)
    loci = pd.read_csv(args.loci, sep="\t", dtype=str, keep_default_na=False)
    variants = load_variants(members)
    if variants["locus_id"].nunique() != len(loci):
        raise ValueError(f"{len(loci)} loci in {args.loci} but {variants['locus_id'].nunique()} in {args.members}")
    log.info("%d loci, %d variants", len(loci), len(variants))

    all_var, all_loops, all_dom, all_loci, qc = [], [], [], [], {"near_kb": args.near_kb, "references": {}}
    for spec in args.reference:
        name, _, d = spec.partition("=")
        d = Path(d)
        files = {k: d / f for k, f in (("loops", "loops.bedpe.gz"), ("domains", "contactDom.bedpe.gz"),
                                        ("compartments", "compartments.bigWig"))}
        missing = [str(p) for p in files.values() if not p.is_file()]
        if missing:
            raise SystemExit(f"reference {name}: missing {missing}")

        loops, loops_meta = read_bedpe(files["loops"])
        domains, dom_meta = read_bedpe(files["domains"])
        lookup, bw_meta = bigwig_lookup(files["compartments"])
        stat_cols = [c for c in loops.columns if re.match(args.loop_stat_regex, c)]
        fdr_cols = [c for c in stat_cols if c.lower().startswith("fdr")]
        if not stat_cols:
            log.warning("%s: no loop columns match %r; header=%s", name, args.loop_stat_regex, loops_meta["columns"])
        log.info("%s: %d loops, %d domains; loop stat columns: %s", name, len(loops), len(domains), stat_cols)

        var, loop_tab, dom_tab = annotate_reference(name, variants, loops, domains, lookup, near_bp, stat_cols)
        all_var.append(var)
        all_loops.append(loop_tab)
        all_dom.append(dom_tab)
        all_loci.append(summarize_loci(loci, var, loop_tab, name, fdr_cols))

        var_chroms = set(variants["chrom"])
        qc["references"][name] = {
            "loops": {k: loops_meta[k] for k in ("path", "n_rows", "columns")},
            "loops_comment_lines": loops_meta["comments"][:10],
            "loops_non_cis": int((loops["_c1"] != loops["_c2"]).sum()),
            "loop_stat_columns": stat_cols,
            "domains": {k: dom_meta[k] for k in ("path", "n_rows", "columns")},
            "domains_non_diagonal": int(((domains["x1"] != domains["y1"]) | (domains["x2"] != domains["y2"])).sum()),
            "compartments": {**bw_meta, "observed_bin_widths": sorted(bw_meta["observed_bin_widths"])},
            "variant_chroms_absent_from_loops": sorted(var_chroms - set(loops["_c1"]) - set(loops["_c2"])),
            "variant_chroms_absent_from_domains": sorted(var_chroms - set(domains["_c1"])),
            "n_variants_without_compartment_score": int(var["compartment_score"].isna().sum()),
            "compartment_label_counts": var["compartment_label"].fillna("missing").value_counts().to_dict(),
            "n_loci_lead_on_anchor": int(all_loci[-1]["lead_on_anchor"].sum()),
            "n_loci_any_member_on_anchor": int((all_loci[-1]["n_members_on_anchor"] > 0).sum()),
        }

    pd.concat(all_var).to_csv(args.out_dir / "variant_annotations.tsv", sep="\t", index=False)
    pd.concat(all_loops).to_csv(args.out_dir / "locus_loops.tsv", sep="\t", index=False)
    pd.concat(all_dom).to_csv(args.out_dir / "locus_domains.tsv", sep="\t", index=False)
    pd.concat(all_loci).to_csv(args.out_dir / "loci_annotated.tsv", sep="\t", index=False)
    with open(args.qc_dir / "annotate_summary.json", "w") as fh:
        json.dump(qc, fh, indent=2, default=str)
    for name, q in qc["references"].items():
        log.info("%s: lead on anchor in %d/%d loci; any member on anchor in %d", name,
                 q["n_loci_lead_on_anchor"], len(loci), q["n_loci_any_member_on_anchor"])


if __name__ == "__main__":
    main()
