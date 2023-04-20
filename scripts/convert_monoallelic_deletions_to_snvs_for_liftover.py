"""
Convert monoallelic deletions to SNVs for liftover
"""

import argparse
import collections
import gzip
import os
import re


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_vcf", help="Path of the input VCF")
    p.add_argument("output_vcf", help="Path of the output VCF")
    args = p.parse_args()

    output_vcf_path = re.sub("[.]b?gz$", "", args.output_vcf)

    counters = collections.defaultdict(int)
    with gzip.open(args.input_vcf, "rt") as f, open(output_vcf_path, "wt") as fo:

        for i, line in enumerate(f):
            if line.startswith("#"):
                fo.write(line)
                continue

            counters["total variants"] += 1
            fields = line.split("\t")
            vcf_ref_allele = fields[3].upper()
            vcf_alt_allele = fields[4].upper()

            # for mono-allelic deletions, convert REF allele to a single nucleotide after saving the original in the INFO field
            if "," not in vcf_alt_allele and len(vcf_ref_allele) > len(vcf_alt_allele) and len(vcf_alt_allele) == 1:
                fields[7] += f";RefAlleleBeforeConversion={vcf_ref_allele}"
                fields[7] += f";AltAlleleBeforeConversion={vcf_alt_allele}"
                fields[3] = fields[3][0]
                fields[4] = "A" if fields[3] != "A" else "C"  # convert to snv
                counters["variants converted to SNVs"] += 1

            fo.write("\t".join(fields))

    os.system(f"bgzip -f {output_vcf_path}")
    os.system(f"tabix -f {output_vcf_path}.gz")

    print(f"Wrote {counters['variants converted to SNVs']} out of {i+1} "
          f"({100*counters['total variants']/(i+1):0.1f}%) lines to {output_vcf_path}.gz")
    print(f"Stats: ")
    for key, count in sorted(counters.items(), key=lambda x: -x[1]):
        print(f"    {count:10,d} ({100*count/counters['total variants']:7.1f}%) {key}")


if __name__ == "__main__":
    main()