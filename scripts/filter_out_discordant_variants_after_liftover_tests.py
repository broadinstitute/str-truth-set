import collections
import pyfaidx
import tempfile
import unittest

from filter_out_discordant_variants_after_liftover import does_one_or_both_alleles_match_t2t_reference_sequence


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        self.temp_fasta_file.write(">chr1\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACACAGACACACTGGGGGGAAAAAAAA\n")
        self.temp_fasta_file.close()
        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def test_does_variant_have_ref_allele(self):
        # fail if VCF ref allele doesn't match the reference
        counter = collections.defaultdict(int)
        true_number_of_AC_repeats = 4
        for pos in 10, 12, 14, 16, 18:
            for num_repeats_offset in -2, -1, 0, 1, 2:
                matches_t2t, num_repeats_in_t2t = does_one_or_both_alleles_match_t2t_reference_sequence(
                    self.fasta_obj, "chr1", pos, ref_base=self.fasta_obj["chr1"][pos], motif="AC",
                    num_repeats_allele1=true_number_of_AC_repeats + num_repeats_offset,
                    num_repeats_allele2=true_number_of_AC_repeats + num_repeats_offset,
                    allow_interruptions=False, repeat_unit_interruption_index=None, counters=counter)
                self.assertTrue(matches_t2t)
                self.assertEqual(num_repeats_in_t2t, true_number_of_AC_repeats)

        for pos in list(range(1, 9)) + [11, 13, 15, 17, 19]:
            matches_t2t, num_repeats_in_t2t = does_one_or_both_alleles_match_t2t_reference_sequence(
                self.fasta_obj, "chr1", pos, ref_base="T", motif="AC",
                num_repeats_allele1=true_number_of_AC_repeats, num_repeats_allele2=true_number_of_AC_repeats + 1,
                allow_interruptions=False, repeat_unit_interruption_index=None, counters=counter)
            self.assertFalse(matches_t2t)
            self.assertEqual(num_repeats_in_t2t, 0)

        matches_t2t, num_repeats_in_t2t = does_one_or_both_alleles_match_t2t_reference_sequence(
            self.fasta_obj, "chr1", 30, ref_base="G", motif="G",
            num_repeats_allele1=6, num_repeats_allele2=6,
            allow_interruptions=False, repeat_unit_interruption_index=None, counters=counter)
        self.assertTrue(matches_t2t)
        self.assertEqual(num_repeats_in_t2t, 6)

    def test_does_variant_have_ref_allele_with_interruptions(self):
        counter = collections.defaultdict(int)
        true_number_of_AC_repeats = 4
        for pos in [10, 14, 18, 22, 26]:
            for num_repeats_offset in -2, -1, 0, 1, 2:
                matches_t2t, num_repeats_in_t2t = does_one_or_both_alleles_match_t2t_reference_sequence(
                    self.fasta_obj, "chr1", pos, ref_base=self.fasta_obj["chr1"][pos], motif="ACAC",
                    num_repeats_allele1=true_number_of_AC_repeats + num_repeats_offset,
                    num_repeats_allele2=true_number_of_AC_repeats + num_repeats_offset,
                    allow_interruptions=True, repeat_unit_interruption_index=1, counters=counter)
                self.assertTrue(matches_t2t)
                self.assertEqual(num_repeats_in_t2t, true_number_of_AC_repeats)
