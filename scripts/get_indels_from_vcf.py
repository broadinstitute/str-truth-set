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
        alt_allele = fields[4]

        if len(ref_allele) != len(alt_allele):
            output_vcf.write(line)

        # update counters
        counters["total variant"] += 1
        if len(ref_allele) != len(alt_allele):
            counters["INDEL"] += 1
            if len(ref_allele) > len(alt_allele):
                counters["DEL"] += 1
            else:
                counters["INS"] += 1
        else:
            counters["SNV/MNV"] += 1


def run(command):
    print(command)
    os.system(command)

run(f"bgzip -f {output_vcf_path}")
run(f"tabix -f {output_vcf_path}.gz")

print(f"Wrote {counters['INDEL']:,d} INDELs to {output_vcf_path}.gz")

for name in "total variant", "SNV/MNV", "INDEL", "INS", "DEL":
    print(f"    {counters[name]:10,d}  {name}s")
