import collections
import pyfaidx
import tempfile
import unittest

from filter_out_discordant_variants_after_liftover import does_variant_have_ref_allele


class Tests(unittest.TestCase):
    def setUp(self):
        self.temp_fasta_file = tempfile.NamedTemporaryFile("w", suffix=".fasta", delete=False)
        self.temp_fasta_file.write(">chr1\n")
        self.temp_fasta_file.write("TGTGTGTGTTACACACACTGTGTGTGTGGGGGG\n")
        self.temp_fasta_file.close()
        self.fasta_obj = pyfaidx.Fasta(self.temp_fasta_file.name, one_based_attributes=False, as_raw=True)

    def test_does_variant_have_ref_allele(self):

        # fail if VCF ref allele doesn't match the reference
        counter = collections.defaultdict(int)
        self.assertFalse(
            does_variant_have_ref_allele(
                self.fasta_obj, "chr1", 1, "A", ["G"], [0,1], counter))
        self.assertDictEqual(counter, {"filtered out: VCF ref allele doesn't match reference sequence": 1})

        # HET (0/1) genotype should always pass regardless of chrom/pos/ref/alt
        counter = collections.defaultdict(int)
        self.assertTrue(
            does_variant_have_ref_allele(
                self.fasta_obj, "chr1", 3, "T", ["TGTGTG"], [0,1], counter))
        self.assertDictEqual(counter, {'kept variants: heterozygous reference genotype': 1})

        # HOM-ALT (1/1) genotype should pass if insertion matches the reference sequence
        counter = collections.defaultdict(int)
        self.assertTrue(
            does_variant_have_ref_allele(
                self.fasta_obj, "chr1", 1, "T", ["TGTGTG"], [1,1], counter))
        self.assertDictEqual(counter, {'kept variants: insertion matches the adjacent reference sequence': 1})

        # multiallelic (1/2) genotype should pass if insertion matches the reference sequence
        counter = collections.defaultdict(int)
        self.assertTrue(
            does_variant_have_ref_allele(
                self.fasta_obj, "chr1", 1, "T", ["TGTGTG", "A"], [1,2], counter))
        self.assertDictEqual(counter, {'kept variants: insertion matches the adjacent reference sequence': 1})

        # HOM-ALT (1/1) genotype should fail if insertion doesn't match the reference sequence
        counter = collections.defaultdict(int)
        self.assertFalse(
            does_variant_have_ref_allele(
                self.fasta_obj, "chr1", 3, "T", ["TATA"], [1,1], counter))
        self.assertDictEqual(counter, {'filtered out variants: INS 1/1': 1})
