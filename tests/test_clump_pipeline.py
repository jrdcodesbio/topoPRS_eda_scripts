"""Tests for the clumping pipeline. All data here is SYNTHETIC.

Unit tests cover the helpers. The end-to-end test (skipped when plink2 is not on PATH) builds
a small phased reference with one perfect-LD block and one independent SNP, then runs
prep_clump_input -> plink2 --clump -> summarize_clumps.

Run from the repo root:  python -m unittest discover -s tests -v
"""
from __future__ import annotations

import json
import random
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

# The repo root is topoEDA/scripts/ on the HPC, so the pipeline scripts live at the repo root.
SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

from gwas_common import allele_pair_key, normalize_chrom, read_pvar_header  # noqa: E402
from summarize_clumps import parse_sp2  # noqa: E402

PLINK2 = shutil.which("plink2")

GWAS_HEADER = ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
               "standard_error", "effect_allele_frequency", "p_value", "rsid",
               "hm_coordinate_conversion", "hm_code", "variant_id"]


class HelperTests(unittest.TestCase):
    def test_normalize_chrom(self):
        s = pd.Series(["chr1", "1", "Chr22", "chrX", "23", "chrM", " 7 "])
        self.assertEqual(normalize_chrom(s).tolist(), ["1", "1", "22", "X", "X", "MT", "7"])

    def test_allele_pair_key_is_order_independent(self):
        k = allele_pair_key(pd.Series(["G", "a", "C"]), pd.Series(["A", "G", "T"]))
        self.assertEqual(k.tolist(), ["A:G", "A:G", "C:T"])

    def test_parse_sp2(self):
        self.assertEqual(parse_sp2("."), [])
        self.assertEqual(parse_sp2("1:200:A:G,1:300:C:T"), ["1:200:A:G", "1:300:C:T"])
        self.assertEqual(parse_sp2("1:200:A:G(1),1:300:C:T(1)"), ["1:200:A:G", "1:300:C:T"])

    def test_read_pvar_header(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.pvar"
            p.write_text("##fileformat=VCFv4.2\n##contig=<ID=1>\n#CHROM\tPOS\tID\tREF\tALT\n1\t10\tv\tA\tG\n")
            cols, skip = read_pvar_header(p)
            self.assertEqual(cols, ["CHROM", "POS", "ID", "REF", "ALT"])
            self.assertEqual(skip, 3)


def write_phased_vcf(path: Path, n_samples: int = 200, seed: int = 1) -> None:
    """Chromosome 1 only. SNPs at 1000/1100/1200 have identical haplotypes (r² = 1).
    The SNP at 50000 is drawn independently. The SNP at 900000 is a second independent locus."""
    rng = random.Random(seed)
    block = [[rng.random() < 0.4 for _ in range(2)] for _ in range(n_samples)]
    indep = [[rng.random() < 0.3 for _ in range(2)] for _ in range(n_samples)]
    far = [[rng.random() < 0.5 for _ in range(2)] for _ in range(n_samples)]

    def gt(h):
        return "|".join("1" if a else "0" for a in h)

    variants = [  # (pos, ref, alt, haplotypes)
        (1000, "A", "G", block),
        (1100, "C", "T", block),
        (1200, "G", "A", block),
        (50000, "T", "C", indep),
        (900000, "A", "C", far),
    ]
    lines = ["##fileformat=VCFv4.2", "##contig=<ID=chr1,length=1000000>",
             "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
             "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
                       + [f"S{i}" for i in range(n_samples)])]
    for pos, ref, alt, haps in variants:
        lines.append("\t".join(["chr1", str(pos), ".", ref, alt, ".", "PASS", ".", "GT"] + [gt(h) for h in haps]))
    path.write_text("\n".join(lines) + "\n")


def write_gwas(path: Path) -> None:
    rows = [
        # Lead of the LD block. Alleles are given as effect/other, opposite to the reference REF/ALT order.
        ["1", "1000", "G", "A", "0.2", "0.03", "0.4", "1e-12", "rs1", "lo", "10", "1_1000_A_G"],
        ["1", "1100", "C", "T", "0.2", "0.03", "0.4", "3e-10", "rs2", "lo", "10", "1_1100_C_T"],
        ["1", "1200", "A", "G", "0.2", "0.03", "0.4", "4e-9", "rs3", "lo", "10", "1_1200_G_A"],
        # Independent of the block, 49 kb away: must form its own clump.
        ["1", "50000", "C", "T", "0.1", "0.02", "0.3", "2e-9", "rs4", "lo", "10", "1_50000_T_C"],
        # Absent from the reference, inside the block window: reported as unmatched, near L0001.
        ["1", "1500", "T", "G", "0.1", "0.02", "0.3", "1e-9", "rs5", "lo", "10", "1_1500_G_T"],
        # Genome-wide significant, but its alleles (G/T) don't match the reference record there (A/C): unmatched.
        ["1", "900000", "G", "T", "0.1", "0.02", "0.5", "1e-8", "rs6", "lo", "10", "1_900000_G_T"],
        # Not significant / missing P: excluded before matching.
        ["1", "900000", "A", "C", "0.01", "0.02", "0.5", "0.4", "rs7", "lo", "10", "1_900000_A_C"],
        ["1", "2000", "A", "C", "0.01", "0.02", "0.5", "NA", "rs8", "lo", "10", "1_2000_A_C"],
    ]
    pd.DataFrame(rows, columns=GWAS_HEADER).to_csv(path, sep="\t", index=False, compression="gzip")


@unittest.skipIf(PLINK2 is None, "plink2 not on PATH")
class EndToEndTest(unittest.TestCase):
    def run_cmd(self, cmd):
        r = subprocess.run([str(c) for c in cmd], capture_output=True, text=True)
        if r.returncode != 0:
            self.fail(f"command failed: {cmd}\nstdout:\n{r.stdout}\nstderr:\n{r.stderr}")
        return r

    def test_pipeline(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            write_phased_vcf(d / "ref.vcf")
            write_gwas(d / "gwas.tsv.gz")
            self.run_cmd([PLINK2, "--vcf", d / "ref.vcf", "--set-all-var-ids", "@:#:$r:$a",
                          "--make-pgen", "--out", d / "ref"])

            py = sys.executable
            self.run_cmd([py, SCRIPTS / "inspect_gwas.py", "--gwas", d / "gwas.tsv.gz", "--out-dir", d / "qc_schema"])
            schema = json.loads((d / "qc_schema" / "summary.json").read_text())
            self.assertEqual(schema["n_rows"], 8)
            self.assertEqual(schema["significance_counts"]["p<5e-08"], 6)
            self.assertEqual(schema["duplicates"]["dup_chrom_pos_rows"], 2)

            self.run_cmd([py, SCRIPTS / "prep_clump_input.py", "--gwas", d / "gwas.tsv.gz",
                          "--ref-pvar", d / "ref.pvar", "--p-max", "5e-8",
                          "--out-clump", d / "clump_input.tsv", "--qc-dir", d / "qc"])
            prep = json.loads((d / "qc" / "prep_summary.json").read_text())
            self.assertEqual(prep["n_sig_complete"], 6)
            self.assertEqual(prep["n_matched_unique_ids"], 4)
            unmatched = pd.read_csv(d / "qc" / "unmatched_sig_snps.tsv", sep="\t", dtype=str)
            self.assertEqual(sorted(unmatched["rsid"]), ["rs5", "rs6"])
            ci = pd.read_csv(d / "clump_input.tsv", sep="\t", dtype=str)
            # plink2's '@' template drops the 'chr' prefix from the input VCF (observed with v2.0.0-a.6.9).
            self.assertEqual(set(ci["ID"]), {"1:1000:A:G", "1:1100:C:T", "1:1200:G:A", "1:50000:T:C"})

            self.run_cmd(["bash", SCRIPTS / "run_clump.sh", PLINK2, d / "ref", d / "clump_input.tsv",
                          d / "out" / "test", "5e-8", "5e-8", "0.1", "250", "1"])
            self.run_cmd([py, SCRIPTS / "summarize_clumps.py", "--clumps", d / "out" / "test.clumps",
                          "--clump-input", d / "clump_input.tsv", "--unmatched", d / "qc" / "unmatched_sig_snps.tsv",
                          "--kb", "250", "--out-dir", d / "res", "--qc-dir", d / "qc"])

            loci = pd.read_csv(d / "res" / "loci.tsv", sep="\t", dtype=str)
            by_lead = loci.set_index("lead_rsid")
            self.assertEqual(sorted(loci["lead_rsid"]), ["rs1", "rs4"])
            self.assertEqual(by_lead.at["rs1", "n_variants"], "3")
            self.assertEqual(by_lead.at["rs1", "span_start"], "1000")
            self.assertEqual(by_lead.at["rs1", "span_end"], "1200")
            self.assertEqual(by_lead.at["rs4", "n_variants"], "1")

            members = pd.read_csv(d / "res" / "locus_members.tsv", sep="\t", dtype=str)
            self.assertEqual(len(members), 4)
            self.assertIn("beta", members.columns)  # original GWAS statistics are kept

            near = pd.read_csv(d / "qc" / "unmatched_near_loci.tsv", sep="\t", dtype=str).set_index("rsid")
            self.assertEqual(near.at["rs5", "within_kb_of_lead"], "True")
            self.assertEqual(near.at["rs5", "nearest_locus"], by_lead.at["rs1", "locus_id"])
            self.assertEqual(near.at["rs6", "within_kb_of_lead"], "False")  # 850 kb from the rs4 lead


if __name__ == "__main__":
    unittest.main()
