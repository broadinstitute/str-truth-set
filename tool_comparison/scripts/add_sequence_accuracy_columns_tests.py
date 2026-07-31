import gzip
import os
import subprocess
import sys
import tempfile
import unittest

import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from add_sequence_accuracy_columns import (
    TOOL_SEQUENCE_COLUMNS, TRUTH_SEQUENCE_COLUMNS,
    compute_sequence_edit_distances, pair_truth_sequences, read_tool_allele_sequences)
from add_concordance_columns import write_alleles_table

SCRIPT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "add_sequence_accuracy_columns.py")

ALLELE_SEQUENCES_HEADER = "\t".join([
    "LocusId", "AlleleSequence: Allele 1", "AlleleSequence: Allele 2",
    "AlleleStatus: Allele 1", "AlleleStatus: Allele 2", "ExtractionStatus"])


class PairTruthSequencesTests(unittest.TestCase):

    def test_both_sequences_available(self):
        # the two sequences are per-haplotype, not size-sorted, so the shorter one becomes Allele 1
        self.assertEqual(pair_truth_sequences("AAAAA", "AAA", 3, 5), ("AAA", "AAAAA"))
        self.assertEqual(pair_truth_sequences("AAA", "AAAAA", 3, 5), ("AAA", "AAAAA"))

    def test_zero_length_allele_is_scored(self):
        # an empty sequence whose reported size is 0 is a real zero-length allele and must be scored, not dropped
        self.assertEqual(pair_truth_sequences("", "AAAAA", 0, 5), ("", "AAAAA"))
        self.assertEqual(pair_truth_sequences("", "", 0, 0), ("", ""))

    def test_missing_sequence_with_nonzero_size_is_unavailable(self):
        # an empty sequence paired with a non-zero reported size means the sequence wasn't stored; scoring it would
        # charge the tool the full length of the allele it actually got right
        self.assertEqual(pair_truth_sequences("", "AAAAA", 3, 5), (None, None))
        # a sequence whose length disagrees with the reported size is not trustworthy either
        self.assertEqual(pair_truth_sequences("AA", "AAAAA", 3, 5), (None, None))

    def test_hemi_locus_duplicates_the_single_sequence(self):
        # HEMI loci (male chrX/chrY outside the PAR) store one sequence plus an empty string, but duplicate the
        # single allele size into both RepeatSize columns. Real example: chrX-2781857-2781866-T, sizes 10/10.
        self.assertEqual(pair_truth_sequences("", "T" * 10, 10, 10), ("T" * 10, "T" * 10))
        self.assertEqual(pair_truth_sequences("T" * 10, "", 10, 10), ("T" * 10, "T" * 10))

        # but a HOM locus that genuinely has a zero-length allele is not treated as HEMI
        self.assertEqual(pair_truth_sequences("", "", 0, 0), ("", ""))


class ComputeSequenceEditDistancesTests(unittest.TestCase):

    def test_index_wise_pairing(self):
        self.assertEqual(compute_sequence_edit_distances("AAA", "AAAAA", "AAA", "AAAAA"), (0, 0))
        # one inserted base on the short allele
        self.assertEqual(compute_sequence_edit_distances("AAA", "AAAAA", "AATA", "AAAAA"), (1, 0))

    def test_unavailable_alleles_score_none(self):
        self.assertEqual(compute_sequence_edit_distances(None, "AAAAA", "AAA", "AAAAA"), (None, 0))
        self.assertEqual(compute_sequence_edit_distances("AAA", "AAAAA", None, None), (None, None))
        # a left-join miss arrives as NaN rather than None, and must be treated the same way
        self.assertEqual(compute_sequence_edit_distances("AAA", "AAAAA", float("nan"), float("nan")), (None, None))

    def test_zero_length_allele_gets_a_real_distance(self):
        # a called zero-length tool allele vs a 5bp truth allele is 5 edits away, not unscorable
        self.assertEqual(compute_sequence_edit_distances("", "AAAAA", "", "AAAAA"), (0, 0))
        self.assertEqual(compute_sequence_edit_distances("AAAAA", "AAAAA", "", "AAAAA"), (5, 0))

    def test_equal_length_tie_break_picks_the_lower_total(self):
        # both truth alleles are 3bp, so either assignment is size-sorted; the pairing with the smaller total wins
        self.assertEqual(compute_sequence_edit_distances("AAA", "CCC", "CCC", "AAA"), (0, 0))
        # ... and the same holds when it's the tool's two alleles that are the same length: index-wise scores
        # 3 + 5 = 8, the swapped assignment scores 0 + 2 = 2
        self.assertEqual(compute_sequence_edit_distances("AAA", "CCCCC", "CCC", "AAA"), (0, 2))

        # an exact tie (index-wise 0 + 4 vs swapped 4 + 0) keeps the index-wise pairing, so the output is deterministic
        self.assertEqual(compute_sequence_edit_distances("AAAA", "AAAA", "AAAA", "CCCC"), (0, 4))

    def test_no_swap_when_neither_side_has_equal_length_alleles(self):
        # both sides are size-sorted with distinct lengths, so the swapped assignment isn't a valid size-sorted
        # pairing even though it would score better (2 + 2 instead of 2 + 4)
        self.assertEqual(compute_sequence_edit_distances("AA", "GGGG", "GG", "AAAA"), (2, 4))


class ReadToolAlleleSequencesTests(unittest.TestCase):

    def write_allele_sequences(self, rows):
        path = os.path.join(self.temp_dir.name, "test.allele_sequences.tsv.gz")
        with gzip.open(path, "wt") as f:
            f.write(ALLELE_SEQUENCES_HEADER + "\n")
            for row in rows:
                f.write("\t".join(row) + "\n")
        return path

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_only_called_alleles_of_ok_loci_are_scored(self):
        path = self.write_allele_sequences([
            ["1-100-110-A", "AAAAA", "AAAAAAA", "called", "called", "ok"],
            # a called zero-length allele (ATaRVa <DEL>) keeps its empty sequence and stays scorable
            ["1-200-210-A", "", "AAAAAAA", "called", "called", "ok"],
            # a no-call allele's empty sequence must not be scored as a zero-length allele
            ["1-300-310-A", "", "", "no_call", "no_call", "ok"],
            ["1-400-410-A", "", "", "no_call", "no_call", "ref_mismatch"],
            ["1-500-510-A", "", "", "no_call", "no_call", "multiallelic"],
        ])
        df = read_tool_allele_sequences(path).set_index("LocusId")

        self.assertEqual(df.loc["1-100-110-A", TOOL_SEQUENCE_COLUMNS[0]], "AAAAA")
        self.assertEqual(df.loc["1-200-210-A", TOOL_SEQUENCE_COLUMNS[0]], "")
        for locus_id in "1-300-310-A", "1-400-410-A", "1-500-510-A":
            for column in TOOL_SEQUENCE_COLUMNS:
                self.assertIsNone(df.loc[locus_id, column], f"{locus_id} {column} should be unscorable")


class MainTests(unittest.TestCase):

    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

        # the truth set table: haplotype-ordered sequences plus the reported per-allele sizes
        self.truth_set_path = os.path.join(self.temp_dir.name, "sample.tandem_repeat_genotypes.tsv.gz")
        pd.DataFrame([
            {"LocusId": "chr1-100-110-A", "Allele1Sequence": "A" * 12, "Allele2Sequence": "A" * 10,
             "RepeatSizeShortAlleleBp": 10, "RepeatSizeLongAlleleBp": 12},
            {"LocusId": "chr1-200-210-A", "Allele1Sequence": "A" * 10, "Allele2Sequence": "A" * 10,
             "RepeatSizeShortAlleleBp": 10, "RepeatSizeLongAlleleBp": 10},
            # a truth locus whose short allele sequence wasn't stored -> both alleles unavailable
            {"LocusId": "chr1-300-310-A", "Allele1Sequence": "", "Allele2Sequence": "A" * 12,
             "RepeatSizeShortAlleleBp": 10, "RepeatSizeLongAlleleBp": 12},
        ]).to_csv(self.truth_set_path, sep="\t", index=False)

        self.allele_sequences_path = os.path.join(self.temp_dir.name, "sample.allele_sequences.tsv.gz")
        with gzip.open(self.allele_sequences_path, "wt") as f:
            f.write(ALLELE_SEQUENCES_HEADER + "\n")
            f.write("\t".join(["1-100-110-A", "A" * 10, "A" * 12, "called", "called", "ok"]) + "\n")
            f.write("\t".join(["1-200-210-A", "A" * 9, "A" * 10, "called", "called", "ok"]) + "\n")
            f.write("\t".join(["1-300-310-A", "A" * 10, "A" * 12, "called", "called", "ok"]) + "\n")

        # the comparison table, as add_tool_results_columns.py writes it (LocusId already chr-stripped)
        self.combined_path = os.path.join(self.temp_dir.name, "combined.with_TRGTv5_results.tsv.gz")
        pd.DataFrame([
            {"LocusId": "1-100-110-A", "MotifSize": 1, "NumRepeats: Allele 1: Truth": 10,
             "NumRepeats: Allele 2: Truth": 12},
            {"LocusId": "1-200-210-A", "MotifSize": 1, "NumRepeats: Allele 1: Truth": 10,
             "NumRepeats: Allele 2: Truth": 10},
            {"LocusId": "1-300-310-A", "MotifSize": 1, "NumRepeats: Allele 1: Truth": 10,
             "NumRepeats: Allele 2: Truth": 12},
            # a locus the tool produced no record for at all
            {"LocusId": "1-400-410-A", "MotifSize": 1, "NumRepeats: Allele 1: Truth": 10,
             "NumRepeats: Allele 2: Truth": 12},
        ]).to_csv(self.combined_path, sep="\t", index=False)

    def tearDown(self):
        self.temp_dir.cleanup()

    def run_script(self, output_path):
        subprocess.run([
            sys.executable, SCRIPT_PATH,
            "--tool", "TRGTv5",
            "--truth-set-genotypes", self.truth_set_path,
            "--allele-sequences", self.allele_sequences_path,
            "--output-tsv", output_path,
            self.combined_path,
        ], check=True, cwd=os.path.dirname(SCRIPT_PATH))
        return pd.read_table(output_path).set_index("LocusId")

    def test_adds_distance_columns_and_drops_every_sequence(self):
        output_path = os.path.join(self.temp_dir.name, "out.tsv.gz")
        df = self.run_script(output_path)

        self.assertEqual(len(df), 4)
        for column in TRUTH_SEQUENCE_COLUMNS + TOOL_SEQUENCE_COLUMNS + [
                "ExtractionStatus", "AlleleStatus: Allele 1", "AlleleStatus: Allele 2"]:
            self.assertNotIn(column, df.columns)
        for column in df.columns:
            self.assertNotIn("Sequence: ", column, "no raw sequence column may survive")

        # an exact match on both alleles
        self.assertEqual(df.loc["1-100-110-A", "SequenceEditDistance: Allele 1: TRGTv5"], 0)
        self.assertEqual(df.loc["1-100-110-A", "SequenceEditDistance: Allele 2: TRGTv5"], 0)
        self.assertEqual(df.loc["1-100-110-A", "SequenceEditDistanceNormalized: Allele 1: TRGTv5"], 0)

        # a 1bp deletion on the short allele: 1 edit over a 10bp truth allele
        self.assertEqual(df.loc["1-200-210-A", "SequenceEditDistance: Allele 1: TRGTv5"], 1)
        self.assertAlmostEqual(df.loc["1-200-210-A", "SequenceEditDistanceNormalized: Allele 1: TRGTv5"], 0.1)
        self.assertEqual(df.loc["1-200-210-A", "SequenceEditDistance: Allele 2: TRGTv5"], 0)

        # the tool called this locus, but the truth sequence is unavailable -> NA, not a distance of 0
        self.assertTrue(pd.isna(df.loc["1-300-310-A", "SequenceEditDistance: Allele 1: TRGTv5"]))
        # the tool produced no record for this locus at all -> NA
        self.assertTrue(pd.isna(df.loc["1-400-410-A", "SequenceEditDistance: Allele 1: TRGTv5"]))

    def test_output_columns_survive_the_alleles_table_melt(self):
        # write_alleles_table() raises a ValueError on any column containing "Allele 1" without "Allele 1: ", so the
        # new columns' names have to be exactly "<Metric>: Allele N: {tool}"
        output_path = os.path.join(self.temp_dir.name, "out.tsv.gz")
        self.run_script(output_path)

        df = pd.read_table(output_path).astype(str)
        alleles_path = os.path.join(self.temp_dir.name, "out.alleles.tsv.gz")
        write_alleles_table(df, alleles_path)

        alleles_df = pd.read_table(alleles_path)
        self.assertEqual(len(alleles_df), 2 * len(df))
        self.assertIn("SequenceEditDistance: Allele: TRGTv5", alleles_df.columns)
        self.assertIn("SequenceEditDistanceNormalized: Allele: TRGTv5", alleles_df.columns)
        # the first locus's two alleles melt into the first and second output rows
        self.assertEqual(alleles_df.iloc[0]["SequenceEditDistance: Allele: TRGTv5"], 0)
        self.assertEqual(alleles_df.iloc[2]["SequenceEditDistance: Allele: TRGTv5"], 1)

    def test_overwrites_the_input_table_by_default(self):
        subprocess.run([
            sys.executable, SCRIPT_PATH,
            "--tool", "TRGTv5",
            "--truth-set-genotypes", self.truth_set_path,
            "--allele-sequences", self.allele_sequences_path,
            self.combined_path,
        ], check=True, cwd=os.path.dirname(SCRIPT_PATH))

        df = pd.read_table(self.combined_path)
        self.assertEqual(len(df), 4)
        self.assertIn("SequenceEditDistance: Allele 1: TRGTv5", df.columns)


if __name__ == "__main__":
    unittest.main()
