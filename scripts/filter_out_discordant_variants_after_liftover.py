"""All true CHM1_CHM13_2 variants should have at least 1 allele that matches the T2T reference after liftover from hg38.
Therefore, the following liftover results point to bad variants in the SynDip Benchmark:
1. variants (insertions or deletions) that have a HOM ALT (ie. homozygous non-reference) genotype after liftover
2. variants that are multi-allelic (eg. have a genotype like "1/2")

Variants that have a heterozgyous genotype (eg. "0/1") are by definition ok since a variant that's encoded as
heterozygous in the VCF has a reference allele that matches the reference.

One additional detail is that, when GATK LiftOver processes an insertion, it doesn't check whether, after liftover,
the insertion matches the new reference (meaning it's now a ref allele). This script performs that check by comparing
the inserted bases to the T2T reference sequence immediately after the insertion.

It also confirms that the ref allele after liftover truely matches the new reference.
"""

import argparse
import collections
import gzip
import os
import pyfaidx
import re

def get_reference_sequence(ref_fasta_obj, chrom, pos_1based, sequence_length):
    """Return the reference sequence string at the given chromosome and position.

    Args:
        ref_fasta_obj (obj): A pyfaidx.Fasta object for reading a reference genome fasta
        chrom (str): chromosome name
        pos_1based (int): 1-based coordinate of the 1st base of the requested sequence
        sequence_length (int): number of base pairs to return
    Return:
        str: The requested reference sequence
    """

    return str(ref_fasta_obj[chrom][pos_1based-1: pos_1based-1 + sequence_length]).upper()


def does_variant_have_ref_allele(
        ref_fasta_obj, chrom, pos, vcf_ref_allele, vcf_alt_alleles, vcf_genotype_indices, counter=None):
    """Returns true if at least one allele in this variant matches the sequence in the reference genome.

    Args:
        ref_fasta_obj (obj): A pyfaidx.Fasta object for reading a reference genome fasta
        chrom (str): The variant's chromosome name
        pos (int): The variant's 1-based genomic position
        vcf_ref_allele (str): The variant's reference allele
        vcf_alt_alleles (list): The list of alt alleles
        vcf_genotype_indices (list): The parsed list of genotype indices (eg. a genotype of "0/1" should be passed
            in as ["0", "1"])
        counter (dict): Optional dictionary of counters for recording stats about kept or discarded variants
    Return:
        bool: True if at least one of the alleles for the given variant matches the genome reference sequence defined by
             ref_fasta_obj.
    """

    try:
        pos = int(pos)
    except ValueError:
        raise ValueError(f"pos arg is not an integer: {pos}")

    if not isinstance(vcf_alt_alleles, list):
        raise ValueError(f"vcf_alt_alleles arg is not a list: {vcf_alt_alleles}")

    if not isinstance(vcf_genotype_indices, list):
        raise ValueError(f"vcf_genotype_indices arg is not a list: {vcf_genotype_indices}")

    # check assumptions about arg values
    if len(vcf_genotype_indices) != 2:
        raise ValueError(f"Unexpected genotype field has {len(vcf_genotype_indices)} genotype indices: {vcf_genotype_indices}")

    for genotype_index in vcf_genotype_indices:
        if genotype_index != ".":
            try:
                genotype_index = int(genotype_index)
            except ValueError:
                raise ValueError(f"Unexpected genotype field has non-numeric genotype indices {vcf_genotype_indices}")
            if genotype_index > len(vcf_alt_alleles):
                raise ValueError(f"Unexpected genotype field has a genotype index that is larger than the number of "
                                 f"alt alleles {vcf_alt_alleles}: {vcf_genotype_indices}")

    # check that the ref allele in the VCF matches the sequence in the reference genome
    fasta_ref_allele = get_reference_sequence(ref_fasta_obj, chrom, pos, len(vcf_ref_allele))
    if fasta_ref_allele != vcf_ref_allele:
        counter["filtered out: VCF ref allele doesn't match reference sequence"] += 1
        return False

    # handle split genotypes like ./1 or 1/.
    if "." in vcf_genotype_indices:
        if counter is not None:
            counter["kept with WARNING: can't validate genotype with '.': " + str(vcf_genotype_indices)] += 1
        return True

    # convert genotype_indices to integers now that '.' genotypes have been dealt with
    vcf_genotype_indices = [int(genotype_index) for genotype_index in vcf_genotype_indices]

    # handle heterozygous genotypes (eg. 0/1)
    if any(genotype_index == 0 for genotype_index in vcf_genotype_indices):
        if counter is not None:
            counter["kept variants: heterozygous reference genotype"] += 1
        return True

    # Homozygous alt. genotypes or multiallelic genotypes are by default considered discordant with the reference (ie.
    # neither allele matches the reference sequence). However, there's a technical issue with gatk LiftoverVCF for
    # alleles that are insertions - even if the inserted sequence exactly matches the new reference, gatk LiftoverVCF
    # doesn't recode the allele as a reference allele. The code below performs this check and, if the inserted
    # sequence does exactly match the reference sequence adjacent to the insertion position, then the variant is
    # rescued and not discarded.
    variant_type = []
    for i, genotype_index in enumerate(set(vcf_genotype_indices)):
        # get the alt allele sequence corresponding to this genotype index
        vcf_alt_allele = vcf_alt_alleles[genotype_index - 1]
        if len(vcf_alt_allele) > len(vcf_ref_allele):
            # get the reference sequence immediately to the right of this variant (assume the variant was left-aligned)
            fasta_alt_allele = get_reference_sequence(ref_fasta_obj, chrom, pos, len(vcf_alt_allele))
            if fasta_alt_allele == vcf_alt_allele:
                if counter is not None:
                    counter["kept variants: insertion matches the adjacent reference sequence"] += 1
                return True

        # record the variant type for logging
        if len(vcf_alt_allele) > len(vcf_ref_allele):
            variant_type.append("INS")
        elif len(vcf_alt_allele) < len(vcf_ref_allele):
            variant_type.append("DEL")
        elif len(vcf_alt_allele) == len(vcf_ref_allele):
            if len(vcf_alt_allele) == 1:
                variant_type.append("SNV")
            else:
                variant_type.append("MNV")

    if counter is not None:
        counter["filtered out variants: " + ",".join(sorted(set(variant_type))) + " " + "/".join(map(str, vcf_genotype_indices))] += 1

    return False


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-R", "--reference-fasta", help="Reference genome fasta path (this is the target reference to "
                                                   "which all variants have been lifted over)", required=True)
    p.add_argument("post_liftover_vcf", help="Path of a VCF produced by gatk LiftoverVCF")
    p.add_argument("output_vcf", nargs="?", help="Optional path of the VCF to which to write all passing variants")
    args = p.parse_args()

    ref_fasta_obj = pyfaidx.Fasta(args.reference_fasta, one_based_attributes=False, as_raw=True)

    if args.output_vcf:
        output_path = re.sub("[.]b?gz$", "", args.output_vcf)
    else:
        output_path = os.path.basename(args.post_liftover_vcf).replace(".gz", "").replace(".vcf", ".passed_liftover_validation.vcf")

    # Example line:
    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	syndip
    # chr1	248752514	.	M	C	30	.	.	GT:AD	1|1:0,2

    counter = collections.defaultdict(int)
    with gzip.open(args.post_liftover_vcf, "rt") as f, open(output_path, "wt") as fo:

        for i, line in enumerate(f):
            if line.startswith("#"):
                fo.write(line)
                continue

            counter['total variants'] += 1

            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])

            vcf_ref_allele = fields[3].upper()
            vcf_alt_alleles = fields[4].upper().split(",")

            vcf_format_field = fields[8]
            vcf_genotype_field = fields[9]
            vcf_format_field_tokens = vcf_format_field.split(":")
            vcf_genotype_field_tokens = vcf_genotype_field.split(":")
            vcf_genotype_field_dict = dict(zip(vcf_format_field_tokens, vcf_genotype_field_tokens))

            # Check assumptions about genotype format
            if "GT" not in vcf_genotype_field_dict:
                raise ValueError(f"Unexpected genotype field in VCF line #{i+1}. GT not found: {line.strip()}")
            vcf_genotype_indices = re.split(r"[/|\\]", vcf_genotype_field_dict["GT"])

            if does_variant_have_ref_allele(
                    ref_fasta_obj, chrom, pos, vcf_ref_allele, vcf_alt_alleles, vcf_genotype_indices, counter):
                fo.write(line)

    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")

    print(f"Wrote {counter['total variants']} out of {i+1} ({100*counter['total variants']/(i+1):0.1f}%) lines to {output_path}.gz")
    print("Stats: ")
    for key, count in sorted(counter.items(), key=lambda x: -x[1]):
        print(f"    {count:6d} ({100*count/counter['total variants']:5.1f}%) {key}")

if __name__ == "__main__":
    main()