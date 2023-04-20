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

from str_analysis.filter_vcf_to_STR_variants import check_if_allele_is_str, get_flanking_reference_sequences
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_allowing_interruptions
from str_analysis.utils.find_repeat_unit import extend_repeat_into_sequence_without_allowing_interruptions


def get_number_of_repeats_in_reference(fasta_obj, chrom, pos, ref_base, allow_interruptions, repeat_unit, repeat_unit_interruption_index):
    """This function is used to determine the number of repeats with the given repeat unit within the reference sequence
    immediately to the left and right of the given position.

    Args:
        fasta_obj (pyfaidx.Fasta): pyfaidx object for the reference fasta file
        chrom (str): chromosome name
        pos (int): 1-based position
        ref_base (str): reference base at the given position
        allow_interruptions (bool): if True, then the repeat unit can be interrupted by a different base
        repeat_unit (str): repeat unit
        repeat_unit_interruption_index (int): index of the repeat unit base that can be interrupted by a different base

    Return:
        int: total number of repeats found in the reference
    """
    left_flanking_reference_sequence, variant_bases, right_flanking_reference_sequence = get_flanking_reference_sequences(
        fasta_obj, chrom, pos, ref_base, ref_base + repeat_unit, num_flanking_bases=2000)

    if allow_interruptions:
        # check for interrupted repeats
        reversed_repeat_unit_interruption_index = (len(repeat_unit) - 1 - repeat_unit_interruption_index) if repeat_unit_interruption_index is not None else None

        num_pure_repeats_left_flank, num_total_repeats_left_flank, reversed_repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1],
            repeat_unit_interruption_index=reversed_repeat_unit_interruption_index)

        if reversed_repeat_unit_interruption_index is not None:
            # reverse the repeat_unit_interruption_index
            repeat_unit_interruption_index = len(repeat_unit) - 1 - reversed_repeat_unit_interruption_index

        num_pure_repeats_right_flank, num_total_repeats_right_flank, repeat_unit_interruption_index = extend_repeat_into_sequence_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence,
            repeat_unit_interruption_index=repeat_unit_interruption_index)

    else:
        num_total_repeats_left_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit[::-1],
            left_flanking_reference_sequence[::-1])

        num_total_repeats_right_flank = extend_repeat_into_sequence_without_allowing_interruptions(
            repeat_unit,
            right_flanking_reference_sequence)

    num_total_repeats = num_total_repeats_left_flank + num_total_repeats_right_flank

    return num_total_repeats


def does_one_or_both_alleles_match_t2t_reference_sequence(
        ref_fasta_obj, chrom, pos, ref_base, motif, num_repeats_allele1, num_repeats_allele2, allow_interruptions,
        repeat_unit_interruption_index, counters=None):
    """Returns true if at least one allele in this variant matches the sequence in the reference genome.

    Args:
        ref_fasta_obj (obj): A pyfaidx.Fasta object for reading a reference genome fasta
        chrom (str): The variant's chromosome name
        pos (int): The variant's 1-based genomic position
        ref_base (str): The variant's reference base
        motif (str): The repeat unit
        num_repeats_allele1 (int): Number of repeats in the first allele
        num_repeats_allele2 (int): Number of repeats in the second allele
        allow_interruptions (bool): If True, then the repeat unit can be interrupted by a different base
        repeat_unit_interruption_index (int): Index of the repeat unit base that can be interrupted by a different base
        counters (dict): Optional dictionary of counters for recording stats about kept or discarded variants

    Return: 2-tuple
        bool: True if at least one of the alleles for the given variant matches the genome reference sequence defined by
             ref_fasta_obj.
        int: Number of repeats found in the reference sequence immediately to the left and right of the variant position.
    """
    num_repeats_in_t2t = get_number_of_repeats_in_reference(
        ref_fasta_obj, chrom, pos, ref_base,
        allow_interruptions=allow_interruptions,
        repeat_unit=motif,
        repeat_unit_interruption_index=repeat_unit_interruption_index)

    diff1 = num_repeats_allele1 - num_repeats_in_t2t
    diff2 = num_repeats_allele2 - num_repeats_in_t2t

    counters["total variants"] += 1
    if abs(diff1) <= 2 or abs(diff2) <= 2:
        counters["passing variants"] += 1

        return True, num_repeats_in_t2t
    else:
        diff = diff1 if abs(diff1) < abs(diff2) else diff2
        if diff < 0:
            counters["T2T alleles bigger than expected"] += 1
        else:
            counters["T2T alleles smaller than expected"] += 1

        return False, num_repeats_in_t2t


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-p", "--log-prefix", default="", help="Optional prefix string to print at the beginning of every line")
    p.add_argument("-R", "--reference-fasta", help="Reference genome fasta path (this is the target reference to "
                                                   "which all variants have been lifted over)", required=True)
    p.add_argument("post_liftover_vcf", help="Path of a VCF produced by gatk LiftoverVCF")
    p.add_argument("output_vcf", help="Path of the VCF where to write all passing variants")
    p.add_argument("failed_vcf", help="Path of the VCF where to write failing variants")
    args = p.parse_args()

    ref_fasta_obj = pyfaidx.Fasta(args.reference_fasta, one_based_attributes=False, as_raw=True)

    output_vcf_path = re.sub("[.]b?gz$", "", args.output_vcf)
    failed_vcf_path = re.sub("[.]b?gz$", "", args.failed_vcf)

    counters = collections.defaultdict(int)
    with gzip.open(args.post_liftover_vcf, "rt") as f, open(output_vcf_path, "wt") as fo, open(failed_vcf_path, "wt") as ff:

        for i, line in enumerate(f):
            if line.startswith("#"):
                fo.write(line)
                continue

            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])

            vcf_ref_allele = fields[3].upper()
            vcf_alt_alleles = fields[4].upper().split(",")
            info_field_dict = {}
            for info_key_value in fields[7].split(";"):
                info_field_tokens = info_key_value.split("=")
                if len(info_field_tokens) > 1:
                    info_field_dict[info_field_tokens[0]] = info_field_tokens[1]
                else:
                    info_field_dict[info_field_tokens[0]] = True

            motif = info_field_dict["Motif"]
            num_repeats_allele1 = int(float(info_field_dict["NumRepeatsShortAllele"]))
            num_repeats_allele2 = int(float(info_field_dict["NumRepeatsLongAllele"]))
            allow_interruptions = "IsPureRepeat" not in info_field_dict
            repeat_unit_interruption_index = int(info_field_dict["MotifInterruptionIndex"]) if allow_interruptions else None
            #num_repeats_in_hg38 = int(float(info_field_dict["NumRepeatsInReference"]))

            found_match, num_repeats_in_t2t = does_one_or_both_alleles_match_t2t_reference_sequence(
                ref_fasta_obj, chrom, pos, vcf_ref_allele[0], motif, num_repeats_allele1, num_repeats_allele2,
                allow_interruptions, repeat_unit_interruption_index, counters)

            fields[7] += f";NumRepeatsInT2T={num_repeats_in_t2t}"
            if found_match:
                fo.write("\t".join(fields))
            else:
                ff.write("\t".join(fields))

            suffix = "PASSED" if found_match else "FAILED"
            for vcf_alt_allele in vcf_alt_alleles:
                if len(vcf_alt_allele) == len(vcf_ref_allele):
                    counters[f"SNV alleles {suffix}"] += 1
                elif len(vcf_alt_allele) > len(vcf_ref_allele):
                    counters[f"INS alleles {suffix}"] += 1
                elif len(vcf_alt_allele) < len(vcf_ref_allele):
                    counters[f"DEL alleles {suffix}"] += 1
                counters[f"total alleles {suffix}"] += 1
                counters[f"total alleles"] += 1

    os.system(f"bgzip -f {output_vcf_path}")
    os.system(f"tabix -f {output_vcf_path}.gz")

    os.system(f"bgzip -f {failed_vcf_path}")
    os.system(f"tabix -f {failed_vcf_path}.gz")

    print(f"Wrote {counters['passing variants']} out of {i+1} ({100*counters['total variants']/(i+1):0.1f}%) lines to {output_vcf_path}.gz")
    print(f"{args.log_prefix} Stats: ")
    for key, count in sorted(counters.items(), key=lambda x: x[0]):
        if ("PASSED" in key or "FAILED" in key) and "alleles" in key:
            continue
        print(f"{args.log_prefix}    {count:7,d} ({100*count/counters['total variants']:5.1f}%) {key}")
    for key, count in sorted(counters.items(), key=lambda x: x[0]):
        if ("PASSED" in key or "FAILED" in key) and "alleles" in key:
            print(f"{args.log_prefix}    {count:7,d} ({100*count/counters['total alleles']:5.1f}%) {key}")


if __name__ == "__main__":
    main()