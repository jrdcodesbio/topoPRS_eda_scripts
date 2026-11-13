"""Synthetic tests for hic_common.read_bedpe and annotate_loci (no HPC data needed).

Coordinates are chosen to pin the 1-based POS vs 0-based half-open BEDPE convention:
a variant at POS = x1 is one base *before* the anchor, POS = x1+1 is its first base,
and POS = x2 is its last base.
"""
from __future__ import annotations

import gzip
import importlib.util
import math
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

from annotate_loci import (  # noqa: E402
    annotate_reference,
    compartment_annotations,
    load_variants,
    merge_intervals,
    summarize_loci,
)
from hic_common import read_bedpe, unique_names  # noqa: E402

LOOP_HEADER = ["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color",
               "observed", "expectedBL", "fdrBL"]
LOOPS = [
    ["chr1", "1000", "2000", "chr1", "50000", "51000", ".", ".", ".", ".", "0,0,0", "30", "5", "0.001"],
    ["chr1", "100000", "105000", "chr1", "200000", "205000", ".", ".", ".", ".", "0,0,0", "12", "4", "0.01"],
]
DOM_HEADER = ["chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score", "strand1", "strand2", "color", "score"]
DOMAINS = [
    ["chr1", "0", "60000", "chr1", "0", "60000", ".", ".", ".", ".", "0,0,0", "0.5"],
    ["chr1", "40000", "55000", "chr1", "40000", "55000", ".", ".", ".", ".", "0,0,0", "0.7"],
]
# (locus, ID, is_lead, POS)
MEMBERS = [
    ("L0001", "v_b", True, 1001),     # first base of anchor 1 of loop0
    ("L0001", "v_a", False, 1000),    # one base before anchor 1 -> near, not overlap
    ("L0001", "v_c", False, 2000),    # last base of anchor 1
    ("L0001", "v_d", False, 50500),   # anchor 2 of loop0; inside both nested domains
    ("L0001", "v_e", False, 30000),   # between loop0 anchors, >10 kb from both
    ("L0002", "v_f", True, 150000),   # between loop1 anchors
    ("L0002", "v_g", False, 203000),  # anchor 2 of loop1
]


def write_bedpe(path: Path, header: list[str], rows: list[list[str]]) -> None:
    with gzip.open(path, "wt") as fh:
        fh.write("#" + "\t".join(header) + "\n")
        fh.write("# juicer_tools version 1.22.01\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")


def members_df() -> pd.DataFrame:
    return pd.DataFrame({
        "locus_id": [m[0] for m in MEMBERS],
        "ID": [m[1] for m in MEMBERS],
        "is_lead": [str(m[2]) for m in MEMBERS],
        "P": ["1e-10"] * len(MEMBERS),
        "REF_CHROM": ["1"] * len(MEMBERS),
        "REF_POS": [str(m[3]) for m in MEMBERS],
        "rsid": [f"rs_{m[1]}" for m in MEMBERS],
    })


def fake_lookup(chrom: str, p0: int):
    """100 kb bins: positive below 100 kb, negative to 200 kb, no data beyond."""
    if chrom != "1" or p0 >= 200_000:
        return float("nan"), None, None
    start = (p0 // 100_000) * 100_000
    return (1.0 if p0 < 100_000 else -0.5), start, start + 100_000


class BedpeTests(unittest.TestCase):
    def test_unique_names(self):
        self.assertEqual(unique_names(["a", "score", "score", "score"]), ["a", "score", "score_2", "score_3"])

    def test_read_bedpe_header_comments_and_duplicates(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "dom.bedpe.gz"
            write_bedpe(p, DOM_HEADER, DOMAINS)
            df, meta = read_bedpe(p)
            self.assertEqual(meta["columns"][-1], "score_2")
            self.assertEqual(meta["comments"], ["# juicer_tools version 1.22.01"])
            self.assertEqual(len(df), 2)
            self.assertEqual(df["x2"].dtype, "int64")
            self.assertEqual(list(df["_c1"]), ["1", "1"])

    def test_read_bedpe_rejects_missing_header(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "bad.bedpe.gz"
            with gzip.open(p, "wt") as fh:
                fh.write("chr1\t1\t2\tchr1\t3\t4\n")
            with self.assertRaises(ValueError):
                read_bedpe(p)


class AnnotateTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        d = Path(cls.tmp.name)
        write_bedpe(d / "loops.bedpe.gz", LOOP_HEADER, LOOPS)
        write_bedpe(d / "dom.bedpe.gz", DOM_HEADER, DOMAINS)
        cls.loops, _ = read_bedpe(d / "loops.bedpe.gz")
        cls.domains, _ = read_bedpe(d / "dom.bedpe.gz")
        cls.variants = load_variants(members_df())
        cls.var, cls.loop_tab, cls.dom_tab = annotate_reference(
            "lung", cls.variants, cls.loops, cls.domains, fake_lookup, 10_000, ["observed", "expectedBL", "fdrBL"])
        cls.v = cls.var.set_index("ID")

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_anchor_overlap_uses_half_open_zero_based(self):
        self.assertEqual(self.v.at["v_a", "n_anchor_loops"], 0)  # POS = x1 is outside
        self.assertEqual(self.v.at["v_a", "n_near_loops"], 1)
        self.assertEqual(self.v.at["v_b", "n_anchor_loops"], 1)  # POS = x1 + 1 is the first base
        self.assertEqual(self.v.at["v_c", "n_anchor_loops"], 1)  # POS = x2 is the last base
        self.assertEqual(self.v.at["v_d", "n_anchor_loops"], 1)

    def test_spanning_is_not_contact(self):
        self.assertEqual(self.v.at["v_e", "n_spanning_loops"], 1)
        self.assertEqual(self.v.at["v_e", "n_anchor_loops"], 0)
        self.assertEqual(self.v.at["v_e", "n_near_loops"], 0)
        self.assertEqual(self.v.at["v_f", "n_spanning_loops"], 1)

    def test_locus_loop_partner_and_stats(self):
        lt = self.loop_tab.set_index(["locus_id", "loop_id"])
        l1 = lt.loc[("L0001", "loop0")]
        self.assertEqual(l1["relation"], "anchor")
        self.assertEqual(l1["anchors_hit"], "1,2")
        self.assertTrue(l1["both_anchors_in_locus"])
        self.assertTrue(l1["lead_on_anchor"])
        self.assertEqual(l1["fdrBL"], "0.001")  # original value carried unchanged
        l2 = lt.loc[("L0002", "loop1")]
        self.assertFalse(l2["lead_on_anchor"])
        self.assertEqual((l2["partner_chrom"], l2["partner_start"], l2["partner_end"]), ("chr1", 100000, 105000))

    def test_domains_nested_and_boundary_distance(self):
        self.assertEqual(self.v.at["v_d", "n_domains"], 2)
        self.assertEqual(self.v.at["v_d", "smallest_domain_span_bp"], 15000)
        self.assertEqual(self.v.at["v_e", "n_domains"], 1)
        self.assertEqual(self.v.at["v_e", "nearest_domain_boundary"], 40000)
        self.assertEqual(self.v.at["v_e", "dist_to_domain_boundary_bp"], 10001)
        self.assertEqual(self.v.at["v_f", "n_domains"], 0)
        self.assertEqual(set(self.dom_tab["locus_id"]), {"L0001"})

    def test_compartments(self):
        self.assertEqual(self.v.at["v_b", "compartment_sign"], "+")
        self.assertEqual(self.v.at["v_f", "compartment_sign"], "-")
        self.assertEqual(self.v.at["v_f", "compartment_bin_start"], 100000)
        self.assertEqual(self.v.at["v_g", "compartment_sign"], "NA")
        self.assertTrue(pd.isna(self.v.at["v_g", "compartment_score"]))
        self.assertEqual(self.v.at["v_b", "compartment_label"], "A")
        self.assertEqual(self.v.at["v_f", "compartment_label"], "B")
        self.assertTrue(pd.isna(self.v.at["v_g", "compartment_label"]))  # no coverage stays missing

    def test_locus_summary(self):
        loci = pd.DataFrame({"locus_id": ["L0001", "L0002"], "lead_id": ["v_b", "v_f"]})
        s = summarize_loci(loci, self.var, self.loop_tab, "lung", ["fdrBL"]).set_index("locus_id")
        self.assertTrue(s.at["L0001", "lead_on_anchor"])
        self.assertEqual(s.at["L0001", "n_members_on_anchor"], 3)
        self.assertEqual(s.at["L0002", "n_partner_regions"], 1)
        self.assertEqual(s.at["L0002", "partner_regions"], "chr1:100000-105000")
        self.assertTrue(math.isclose(s.at["L0002", "min_fdrBL_anchor_loops"], 0.01))
        self.assertFalse(s.at["L0002", "lead_on_anchor"])
        self.assertEqual(s.at["L0001", "lead_compartment_label"], "A")
        self.assertEqual(s.at["L0001", "n_members_A"], 5)
        self.assertEqual(s.at["L0002", "lead_compartment_label"], "B")
        self.assertEqual((s.at["L0002", "n_members_B"], s.at["L0002", "n_members_compartment_na"]), (1, 1))

    def test_near_window_zero_disables_near(self):
        var, _, _ = annotate_reference("lung", self.variants, self.loops, self.domains, fake_lookup, 0, [])
        self.assertEqual(var.set_index("ID").at["v_a", "n_near_loops"], 0)


class HelperTests(unittest.TestCase):
    def test_merge_intervals(self):
        merged = merge_intervals([("chr1", 10, 20), ("chr1", 15, 30), ("chr1", 30, 40), ("chr2", 0, 5)])
        self.assertEqual(merged, [("chr1", 10, 40), ("chr2", 0, 5)])

    def test_compartment_zero_score(self):
        v = load_variants(members_df().iloc[:1])
        out = compartment_annotations(v, lambda c, p: (0.0, 0, 10))
        self.assertEqual(out["compartment_sign"].iloc[0], "0")
        self.assertEqual(out["compartment_label"].iloc[0], "unassigned")


@unittest.skipIf(importlib.util.find_spec("pyBigWig") is None, "pyBigWig not installed")
class BigWigTests(unittest.TestCase):
    def test_bigwig_lookup_roundtrip(self):
        import pyBigWig

        from annotate_loci import bigwig_lookup

        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "c.bw"
            bw = pyBigWig.open(str(p), "w")
            bw.addHeader([("chr1", 300_000)])
            bw.addEntries(["chr1", "chr1"], [0, 100_000], ends=[100_000, 200_000], values=[1.5, -2.0])
            bw.close()
            lookup, meta = bigwig_lookup(p)
            self.assertEqual(lookup("1", 99_999), (1.5, 0, 100_000))
            self.assertEqual(lookup("1", 100_000), (-2.0, 100_000, 200_000))
            self.assertTrue(math.isnan(lookup("1", 250_000)[0]))
            self.assertTrue(math.isnan(lookup("2", 10)[0]))
            self.assertEqual(meta["observed_bin_widths"], {100_000})


if __name__ == "__main__":
    unittest.main()
