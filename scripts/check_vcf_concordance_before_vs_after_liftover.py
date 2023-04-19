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

p = argparse.ArgumentParser()
p.add_argument("-p", "--log-prefix", default="", help="Optional prefix string to print at the beginning of log lines")
p.add_argument("-o", "--output-vcf", help="Optional output vcf path for variants that pass checks")
p.add_argument("vcf_before_liftover", help="Variant calls before hg38 => T2T => hg38 liftover")
p.add_argument("vcf_after_liftover", help="Variants calls after hg38 => T2T => hg38 liftover")
args = p.parse_args()

if args.output_vcf:
    output_path = re.sub("[.]b?gz$", "", args.output_vcf)
else:
    output_path = os.path.basename(args.vcf_after_liftover).replace(".gz", "").replace(
        ".vcf", ".passed_liftover_checks.vcf")

# Example line:
# #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	syndip
# chr1	248752514	.	M	C	30	.	.	GT:AD	1|1:0,2

counters = collections.defaultdict(int)
alleles_before_liftover = set()
fopen = gzip.open if args.vcf_before_liftover.endswith("gz") else open
with fopen(args.vcf_before_liftover, "rt") as vcf_file_before_liftover:
    for line in vcf_file_before_liftover:
        if line.startswith("#"):
            continue

        fields = line.strip().split("\t")
        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3].upper()
        alt_alleles = [a.upper() for a in fields[4].split(",")]
        counters['TOTAL variants before liftover'] += 1
        for alt_allele in alt_alleles:
            alleles_before_liftover.add((chrom, pos, ref, alt_allele))
            counters['TOTAL alleles before liftover'] += 1

print(f"Parsed {counters['TOTAL variants before liftover']:,d} variants from {args.vcf_before_liftover}")

fopen = gzip.open if args.vcf_after_liftover.endswith("gz") else open
with fopen(args.vcf_after_liftover, "rt") as vcf_file_after_liftover, open(output_path, "wt") as fo:

    for i, line in enumerate(vcf_file_after_liftover):
        if line.startswith("#"):
            fo.write(line)
            continue

        counters["TOTAL variants after liftover"] += 1
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

        # If this variant is a deletion converted to a SNV by the convert_monoallelic_deletinos_to_snvs_for_liftover.py
        # script, prior to hg38 => T2T liftover, then restore the original ref allele.
        if "RefAlleleBeforeConversion" in info_field_dict:
            counters["variants converted from SNV back to mono-allelic DEL"] += 1
            fields[3] = ref = info_field_dict["RefAlleleBeforeConversion"]
            fields[4] = info_field_dict["AltAlleleBeforeConversion"]
            alt_alleles = [fields[4]]

        failed_liftover_due_to_straddles_multiple_intervals_error = fields[6] == "IndelStraddlesMultipleIntevals"
        fields[7] += f";SkippedValidation={failed_liftover_due_to_straddles_multiple_intervals_error}"
        fields[6] = "PASS"      # reset the filter field

        allele_not_found_in_vcf_before_liftover = False
        for alt_allele in alt_alleles:
            counters["TOTAL alleles after liftover"] += 1
            if (chrom, pos, ref, alt_allele) not in alleles_before_liftover:
                allele_not_found_in_vcf_before_liftover = True
                counters["alleles failed"] += 1
                if len(ref) > len(alt_allele):
                    counters["DEL alleles failed"] += 1
                else:
                    counters["INS alleles failed"] += 1

        if allele_not_found_in_vcf_before_liftover:
            counters["variants had a different position after hg38 => T2T => hg38"] += 1
            continue

        counters["variants passed"] += 1
        fo.write("\t".join(fields))

os.system(f"bgzip -f {output_path}")
os.system(f"tabix -f {output_path}.gz")

total_after_liftover = counters['TOTAL variants after liftover']
print(f"Parsed {total_after_liftover:,d} variants from {args.vcf_after_liftover}")

print(f"Wrote {counters['variants passed']:,d} out of {counters['TOTAL variants before liftover']:,d} "
      f"({100*counters['variants passed']/counters['TOTAL variants before liftover']:0.1f}%) "
      f"variants to {output_path}.gz")

for key in "variants had a different position after hg38 => T2T => hg38", "variants converted from SNV back to mono-allelic DEL":
    print(args.log_prefix, f"  {counters[key]:7,d} out of {total_after_liftover:7,d} ({counters[key]/total_after_liftover:6,.1%}) {key}")
