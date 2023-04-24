"""This script converts a VCF file to a TSV file.
It uses pandas to figure out the TSV header and so is intended for small VCFs that easily fit into memory.
"""

import argparse
import gzip
import pandas as pd
import re


def parse_vcf_line(line):
    fields = line.rstrip().split("\t")
    info_field = fields[7]

    info_field_dict = {}
    for info_key_value in info_field.split(";"):
        info_field_tokens = info_key_value.split("=")
        if len(info_field_tokens) > 1:
            info_field_dict[info_field_tokens[0]] = info_field_tokens[1]
        else:
            info_field_dict[info_field_tokens[0]] = True

    return info_field_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("vcf_path", help="VCF input file")
    args = parser.parse_args()

    filename_prefix = re.sub("[.]vcf(.b?gz)?$", "", args.vcf_path)
    fopen = gzip.open if args.vcf_path.endswith("gz") else open

    with fopen(args.vcf_path, "rt") as vcf_file:
        output_rows = [parse_vcf_line(line) for line in vcf_file if not line.startswith("#")]

    df = pd.DataFrame(output_rows)
    df.to_csv(f"{filename_prefix}.variants.tsv.gz", sep="\t", index=False, header=True)

    print(f"Wrote {len(df):,d} rows to {filename_prefix}.variants.tsv.gz with columns:\n    {', '.join(df.columns)}")


if __name__ == "__main__":
    main()