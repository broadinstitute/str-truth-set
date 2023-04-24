"""This script checks whether any variants changed positions or allele sequences due to hg38=>T2T=>hg38 liftover.

It takes 2 VCFs:
1. the original VCF from before hg38=>T2T=>hg38 liftover
2. the VCF with all variants that survived hg38=>T2T=>hg38 liftover

It then checks that all variants in the "after" VCF are also present in the "before" VCF. If a variant in the
"after" VCF is new relative to the "before" VCF, this means the position or alleles changed due to hg38=>T2T=>hg38.
"""

import argparse
import collections
import gzip
import os
import re

COMMON_TSV_HEADER = [
    "LocusId",
    "Locus",
    "Motif",
    "NumRepeatsInReference",
    "IsPureRepeat",
]


VARIANTS_TSV_HEADER = COMMON_TSV_HEADER + ["NumRepeatsShortAllele", "NumRepeatsLongAllele"]
ALLELES_TSV_HEADER = COMMON_TSV_HEADER + ["NumRepeats"]


def parse_vcf_line(line):

    fields = line.split("\t")
    chrom = fields[0]
    pos = int(fields[1])
    ref = fields[3].upper()
    alt_alleles = [a.upper() for a in fields[4].split(",")]

    # Parse the INFO field
    info_field_dict = {}
    for info_key_value in fields[7].split(";"):
        info_field_tokens = info_key_value.split("=")
        if len(info_field_tokens) > 1:
            info_field_dict[info_field_tokens[0]] = info_field_tokens[1]
        else:
            info_field_dict[info_field_tokens[0]] = True

    return fields, chrom, pos, ref, alt_alleles, info_field_dict


def save_variants_that_passed_validation(args, alleles_before_validation, counters):
    if args.output_vcf:
        output_path = re.sub("[.]b?gz$", "", args.output_vcf)
    else:
        output_path = os.path.basename(args.vcf_after_validation).replace(".gz", "").replace(
            ".vcf", ".passed_validation_checks.vcf")

    alleles_passed_validation = set()
    fopen = gzip.open if args.vcf_after_validation.endswith("gz") else open
    with fopen(args.vcf_after_validation, "rt") as vcf_file_after_validation, open(output_path, "wt") as fo:

        for i, line in enumerate(vcf_file_after_validation):
            if line.startswith("#"):
                fo.write(line)
                continue

            counters["TOTAL variants after validation"] += 1
            fields, chrom, pos, ref, alt_alleles, info_field_dict = parse_vcf_line(line)

            # If this variant is a deletion converted to a SNV by the convert_monoallelic_deletinos_to_snvs_for_validation.py
            # script, prior to hg38 => T2T validation, then restore the original ref allele.
            if "RefAlleleBeforeConversion" in info_field_dict:
                counters["variants converted from SNV back to mono-allelic DEL"] += 1
                fields[3] = ref = info_field_dict["RefAlleleBeforeConversion"]
                fields[4] = info_field_dict["AltAlleleBeforeConversion"]
                alt_alleles = [fields[4]]

            failed_validation_due_to_straddles_multiple_intervals_error = fields[6] == "IndelStraddlesMultipleIntevals"
            fields[7] += f";SkippedValidation={failed_validation_due_to_straddles_multiple_intervals_error}"
            fields[6] = "PASS"      # reset the filter field

            allele_not_found_in_vcf_before_validation = False
            for alt_allele in alt_alleles:
                counters["TOTAL alleles after validation"] += 1
                if (chrom, pos, ref, alt_allele) not in alleles_before_validation:
                    allele_not_found_in_vcf_before_validation = True
                    counters["alleles failed"] += 1
                    if len(ref) > len(alt_allele):
                        counters["DEL alleles failed"] += 1
                    else:
                        counters["INS alleles failed"] += 1

            if allele_not_found_in_vcf_before_validation:
                counters["variants had a different position after hg38 => T2T => hg38"] += 1
                continue

            if args.failed_validation_output_prefix:
                # record all alleles that passed validation
                for alt_allele in alt_alleles:
                    alleles_passed_validation.add((chrom, pos, ref, alt_allele))

            counters["variants passed"] += 1
            fo.write("\t".join(fields))

    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")

    total_variants_after_validation = counters['TOTAL variants after validation']
    print(args.log_prefix, f"Parsed {total_variants_after_validation:,d} variants from {args.vcf_after_validation}")

    print(args.log_prefix, f"Wrote {counters['variants passed']:,d} out of {counters['TOTAL variants before validation']:,d} "
          f"({100*counters['variants passed']/counters['TOTAL variants before validation']:0.1f}%) "
          f"variants to {output_path}.gz")

    for key in "variants had a different position after hg38 => T2T => hg38", "variants converted from SNV back to mono-allelic DEL":
        print(args.log_prefix, f"  {counters[key]:7,d} out of {total_variants_after_validation:7,d} "
                               f"({counters[key]/total_variants_after_validation:0,.1%}) {key}")

    total_variants_before_validation = counters['TOTAL variants before validation']
    total_alleles_before_validation = counters['TOTAL alleles before validation']
    total_alleles_after_validation = counters['TOTAL alleles after validation']
    print(args.log_prefix, f"Found {len(alleles_passed_validation):,d} out of {total_alleles_before_validation:,d} ("
          f"{len(alleles_passed_validation)/total_alleles_before_validation:0,.1%}) "
          f"alleles passed validation, while "
          f"{total_variants_before_validation - total_variants_after_validation:,d} ("
          f"{(total_variants_before_validation - total_variants_after_validation)/total_variants_before_validation:0.1%}) "
          f"variants and "
          f"{total_alleles_before_validation - total_alleles_after_validation:,d} ("
          f"{(total_alleles_before_validation - total_alleles_after_validation)/total_alleles_before_validation:0.1%}) "
          f"alleles failed validation.")
    return alleles_passed_validation


def save_variants_that_failed_validation(args, alleles_passed_validation, counters):
    fopen = gzip.open if args.vcf_before_validation.endswith("gz") else open
    with fopen(args.vcf_before_validation, "rt") as vcf_file_before_validation, \
            open(f"{args.failed_validation_output_prefix}.vcf", "wt") as fo_vcf, \
            open(f"{args.failed_validation_output_prefix}.variants.tsv", "wt") as fo_variants_tsv, \
            open(f"{args.failed_validation_output_prefix}.alleles.tsv", "wt") as fo_alleles_tsv:

        fo_variants_tsv.write("\t".join(VARIANTS_TSV_HEADER) + "\n")
        fo_alleles_tsv.write("\t".join(ALLELES_TSV_HEADER) + "\n")

        for i, line in enumerate(vcf_file_before_validation):
            if line.startswith("#"):
                fo_vcf.write(line)
                continue

            fields, chrom, pos, ref, alt_alleles, info_field_dict = parse_vcf_line(line)
            if all((chrom, pos, ref, alt_allele) in alleles_passed_validation for alt_allele in alt_alleles):
                # this variant passed validation
                continue

            counters["failed variants"] += 1
            fo_vcf.write(line)


            # Write INFO to TSV
            fo_variants_tsv.write("\t".join([info_field_dict[key] for key in VARIANTS_TSV_HEADER]) + "\n")

            for alt_allele_key in "NumRepeatsShortAllele", "NumRepeatsLongAllele":
                counters["failed alleles"] += 1
                info_field_dict["NumRepeats"] = info_field_dict[alt_allele_key]
                fo_alleles_tsv.write("\t".join([info_field_dict[key] for key in ALLELES_TSV_HEADER]) + "\n")

    os.system(f"bgzip -f {args.failed_validation_output_prefix}.variants.tsv")
    os.system(f"bgzip -f {args.failed_validation_output_prefix}.alleles.tsv")
    os.system(f"bgzip -f {args.failed_validation_output_prefix}.vcf")
    os.system(f"tabix -f {args.failed_validation_output_prefix}.vcf.gz")

    print(args.log_prefix, "Wrote", f"{counters['failed variants']:,d}", "variants that failed validation to",
          f"{args.failed_validation_output_prefix}.vcf.gz and {args.failed_validation_output_prefix}.variants.tsv.gz")
    print(args.log_prefix, "Wrote", f"{counters['failed alleles']:,d}", "alleles that failed validation to",
          f"{args.failed_validation_output_prefix}.alleles.tsv.gz")


def main():

    p = argparse.ArgumentParser()
    p.add_argument("-p", "--log-prefix", default="", help="Optional prefix string to print at the beginning of log lines")
    p.add_argument("-o", "--output-vcf", help="Optional output vcf path for variants that pass checks")
    p.add_argument("-f", "--failed-validation-output-prefix", help="Optional output path prefix for all variants that "
                                                                   "failed validation")
    p.add_argument("vcf_before_validation", help="Variant calls before hg38 => T2T => hg38 liftover")
    p.add_argument("vcf_after_validation", help="Variants calls after hg38 => T2T => hg38 liftover")
    args = p.parse_args()

    # Example line:
    # #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	syndip
    # chr1	248752514	.	M	C	30	.	.	GT:AD	1|1:0,2

    counters = collections.defaultdict(int)
    alleles_before_validation = set()
    fopen = gzip.open if args.vcf_before_validation.endswith("gz") else open
    with fopen(args.vcf_before_validation, "rt") as vcf_file_before_validation:
        for line in vcf_file_before_validation:
            if line.startswith("#"):
                continue
            counters['TOTAL variants before validation'] += 1
            fields, chrom, pos, ref, alt_alleles, _ = parse_vcf_line(line)
            for alt_allele in alt_alleles:
                alleles_before_validation.add((chrom, pos, ref, alt_allele))
                counters['TOTAL alleles before validation'] += 1

    print(args.log_prefix, f"Parsed {counters['TOTAL variants before validation']:,d} variants from {args.vcf_before_validation}")
    print(args.log_prefix, f"Parsed {counters['TOTAL alleles before validation']:,d} alleles from {args.vcf_before_validation}")

    alleles_that_passed_validation = save_variants_that_passed_validation(args, alleles_before_validation, counters)

    if args.failed_validation_output_prefix:
        save_variants_that_failed_validation(args, alleles_that_passed_validation, counters)

    print(args.log_prefix, "Done")


if __name__ == "__main__":
    main()