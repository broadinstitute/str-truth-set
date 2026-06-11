"""Tests for compute_truth_set_tsv_for_comparisons.py, focused on normalizing a v2 tandem_repeat_genotypes tsv
(from str_analysis.filter_vcf_to_tandem_repeats genotype) into the for_comparison table used by the accuracy plots.
"""

import gzip
import os
import subprocess
import tempfile
import unittest

import pandas as pd

SCRIPT_PATH = os.path.join(os.path.dirname(__file__), "compute_truth_set_tsv_for_comparisons.py")

# Minimal subset of the genotype tsv columns (the rest are filled with empty strings below).
GENOTYPE_COLUMNS = [
    "Chrom", "Start0Based", "End", "Locus", "LocusId", "Motif", "CanonicalMotif", "MotifSize",
    "NumRepeatsInReference", "NumRepeatsShortAllele", "NumRepeatsLongAllele",
    "RepeatSizeShortAlleleBp", "RepeatSizeLongAlleleBp", "Zygosity", "IsPureRepeat", "RepeatPurity",
    "RepeatPurityShortAllele", "RepeatPurityLongAllele", "NumOverlappingVariants", "VariantPositions",
]

# MotifSize 3 (CAG), reference has 4 repeats (12 bp) spanning chr1:100-112 (0-based start).
# Each row: (locus_suffix, num_short, num_long, zygosity, purity_short, purity_long)
GENOTYPE_ROWS = [
    ("het",     4, 5, "HET",  1.0,  0.9),    # ref/alt -> HET; short allele IS the reference allele
    ("homalt",  6, 6, "HOM",  0.95, 0.95),   # homozygous alt -> HOM
    ("multi",   3, 6, "HET",  0.8,  0.7),    # two alt alleles, neither == ref -> MULTI
    ("homref",  4, 4, "HOM",  1.0,  1.0),    # non-variant (both == ref) -> dropped from truth set
    ("hemi",    5, 5, "HEMI", 0.85, 0.85),   # single allele -> HEMI
]


def _build_genotype_tsv(path):
    rows = []
    for suffix, num_short, num_long, zygosity, purity_short, purity_long in GENOTYPE_ROWS:
        row = {c: "" for c in GENOTYPE_COLUMNS}
        row.update({
            "Chrom": "chr1",
            "Start0Based": 100,
            "End": 112,
            "Locus": f"chr1:101-112",
            "LocusId": f"chr1-100-112-CAG-{suffix}",
            "Motif": "CAG",
            "CanonicalMotif": "ACG",
            "MotifSize": 3,
            "NumRepeatsInReference": 4,
            "NumRepeatsShortAllele": num_short,
            "NumRepeatsLongAllele": num_long,
            "RepeatSizeShortAlleleBp": num_short * 3,
            "RepeatSizeLongAlleleBp": num_long * 3,
            "Zygosity": zygosity,
            "IsPureRepeat": purity_short > 0.99 and purity_long > 0.99,
            "RepeatPurity": min(purity_short, purity_long),
            "RepeatPurityShortAllele": purity_short,
            "RepeatPurityLongAllele": purity_long,
            "NumOverlappingVariants": 1,
            "VariantPositions": "101",
        })
        rows.append(row)
    pd.DataFrame(rows, columns=GENOTYPE_COLUMNS).to_csv(path, sep="\t", index=False)


class TestNormalizeV2GenotypeFormat(unittest.TestCase):

    def test_for_comparison_from_v2_genotype_tsv(self):
        with tempfile.TemporaryDirectory() as tmp_dir:
            genotype_tsv = os.path.join(tmp_dir, "SAMPLE1.tandem_repeat_genotypes.tsv")
            _build_genotype_tsv(genotype_tsv)

            subprocess.run(
                ["python3", SCRIPT_PATH, "--output-dir", tmp_dir, genotype_tsv],
                check=True, capture_output=True, text=True)

            out_path = os.path.join(tmp_dir, "SAMPLE1.tandem_repeat_genotypes.for_comparison.tsv")
            self.assertTrue(os.path.exists(out_path), f"missing output: {out_path}")
            df = pd.read_table(out_path)
            df = df.set_index(df["LocusId"].str.rsplit("-", n=1).str[-1])  # index by the locus suffix

            # the non-variant (hom-ref) locus is dropped, leaving the 4 variant loci
            self.assertEqual(set(df.index), {"het", "homalt", "multi", "hemi"})

            # genotype classification relative to the reference allele
            self.assertEqual(df.loc["het", "HET_or_HOM_or_HEMI_or_MULTI"], "HET")
            self.assertEqual(df.loc["homalt", "HET_or_HOM_or_HEMI_or_MULTI"], "HOM")
            self.assertEqual(df.loc["multi", "HET_or_HOM_or_HEMI_or_MULTI"], "MULTI")
            self.assertEqual(df.loc["hemi", "HET_or_HOM_or_HEMI_or_MULTI"], "HEMI")
            self.assertTrue(df.loc["multi", "IsMultiallelic"])
            self.assertFalse(df.loc["het", "IsMultiallelic"])

            # the key fix: the reference allele (Allele 1 of the het ref/alt locus) carries a real repeat purity
            self.assertAlmostEqual(df.loc["het", "RepeatPurity: Allele 1"], 1.0)
            self.assertAlmostEqual(df.loc["het", "RepeatPurity: Allele 2"], 0.9)
            # no NaN/blank purity anywhere
            self.assertFalse(df["RepeatPurity: Allele 1"].isna().any())
            self.assertFalse(df["RepeatPurity: Allele 2"].isna().any())


if __name__ == "__main__":
    unittest.main()
