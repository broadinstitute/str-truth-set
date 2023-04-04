"""This script takes a vcf path and writes outs a vcf that contains just the INDELs"""

import argparse
import collections
import gzip
import os
import re

p = argparse.ArgumentParser()
p.add_argument("--output-vcf", help="Path of output vcf")
p.add_argument("vcf_path")
args = p.parse_args()

output_vcf_path = args.output_vcf or (
        re.sub(".vcf(.b?gz)?$", "", args.vcf_path) + ".INDELs.vcf")

fopen = gzip.open if args.vcf_path.endswith("gz") else open
counters = collections.defaultdict(int)
print(f"Parsing {args.vcf_path}")
with fopen(args.vcf_path, "rt") as input_vcf, open(output_vcf_path, "wt") as output_vcf:
    for line in input_vcf:
        if line.startswith("#"):
            output_vcf.write(line)
            continue

        fields = line.split("\t")
        ref_allele = fields[3]
        alt_alleles = fields[4].split(",")

        has_indel_allele = False
        for alt_allele in alt_alleles:
            if alt_allele == "*":
                continue

            counters["total allele"] += 1
            if len(ref_allele) == len(alt_allele):
                counters["SNV/MNV"] += 1
            else:
                if len(ref_allele) < len(alt_allele) and alt_allele.startswith(ref_allele):
                    has_indel_allele = True
                    counters["INS allele"] += 1
                elif len(ref_allele) > len(alt_allele) and ref_allele.startswith(alt_allele):
                    has_indel_allele = True
                    counters["DEL allele"] += 1

        counters["total variant"] += 1
        if has_indel_allele:
            counters["INDEL"] += 1
            output_vcf.write(line)


def run(command):
    print(command)
    os.system(command)

run(f"bgzip -f {output_vcf_path}")
run(f"tabix -f {output_vcf_path}.gz")

print(f"Wrote {counters['INDEL']:,d} INDELs to {output_vcf_path}.gz")

# print stats
for name in "total variant", "SNV/MNV", "INDEL":
    print(f"  {counters[name]:10,d} out of {counters['total variant']:10,d} "
          f"({100*counters[name]/counters['total variant']:6.1f}%) {name} variants")

for name in "total allele", "INDEL allele", "INS allele", "DEL allele":
    if name == "INDEL allele":
        value = counters["INS allele"] + counters["DEL allele"]
    else:
        value = counters[name]

    print(f"  {value:10,d} out of {counters['total allele']:10,d} ({100*value/counters['total allele']:6.1f}%) {name}s")
