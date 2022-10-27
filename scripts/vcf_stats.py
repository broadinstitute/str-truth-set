"""This script takes a VCF file path and prints a summary of the number of variants and alleles"""

import argparse
import collections
import gzip
import sys

p = argparse.ArgumentParser()
p.add_argument("-m", "--min-percent", default=0, type=float,
               help="Only print stats (eg. % of variants that is multi-allelic) if the % is greater than this threshold")
p.add_argument("-l", "--label", help="Optional text label to describe the VCF")
p.add_argument("vcf_path")
args = p.parse_args()

if args.label:
    print(f"VCF stats for {args.label}",  "-- "*50, args.vcf_path)
else:
    print(f"VCF stats for ",  "-- "*50, args.vcf_path)

# step 1: Compute stats
counters = collections.defaultdict(int)
fopen = gzip.open if args.vcf_path.endswith("gz") else open
for i, line in enumerate(fopen(args.vcf_path, "rt")):
    if line.startswith("#"):
        continue

    # VCF Header:  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  syndip
    fields = line.split("\t")
    try:
        chrom = fields[0]
        pos = fields[1]
        ref = fields[3]
        alt_alleles = fields[4].split(",")
        fltr = fields[6]
        info = fields[7]
        genotype_fields = fields[9].split(":")
    except Exception as e:
        raise ValueError(f"Unable to parse line #{i+1}: '{line}' in {args.vcf_path}. {e}")

    if len(alt_alleles) == 0:
        raise ValueError(f"No alt alleles found in line #{i+1}: '{line}' in {args.vcf_path}")

    counters["TOTAL variants"] += 1
    counters[f"genotype: {genotype_fields[0]}"] += 1

    # count alleles
    insertion_counter = 0
    deletion_counter = 0
    snv_counter = 0
    for alt in alt_alleles:
        counters["TOTAL alleles"] += 1

        if len(ref) == len(alt):
            snv_counter += 1
            counters[f"SNV alleles"] += 1
            if len(ref) > 1:
                counters[f"SNV alleles (MNV)"] += 1
        elif len(ref) > len(alt):
            deletion_counter += 1
            counters[f"DEL alleles"] += 1
        elif len(ref) < len(alt):
            insertion_counter += 1
            counters[f"INS alleles"] += 1
        else:
            raise ValueError("Invalid state")

    # count variant types
    if insertion_counter > 0 and deletion_counter == 0 and snv_counter == 0:
        counters["INS variants"] += 1
        if insertion_counter > 1: counters["INS variants (multiallelic)"] += 1
    elif insertion_counter == 0 and deletion_counter > 0 and snv_counter == 0:
        counters["DEL variants"] += 1
        if deletion_counter > 1: counters["DEL variants (multiallelic)"] += 1
    elif insertion_counter == 0 and deletion_counter == 0 and snv_counter > 0:
        counters["SNV variants"] += 1
        if snv_counter > 1: counters["SNV variants (multiallelic)"] += 1
    else:
        var_types = []
        if insertion_counter > 0: var_types.append("INS")
        if deletion_counter > 0:  var_types.append("DEL")
        if snv_counter > 0:       var_types.append("SNV")
        var_types = "/".join(var_types)
        counters[f"mixed {var_types} variants (multiallelic)"] += 1


# step 2: Print totals
print("{:10,d}".format(counters["TOTAL variants"]), "TOTAL variants")
print(
    "{:10,d}".format(counters["TOTAL alleles"]), "TOTAL alleles",
    ("(%0.2f alleles per variant)" % (counters["TOTAL alleles"]/counters["TOTAL variants"])) if counters["TOTAL variants"] > 0 else ""
)

# step 3: Print all counters
for key, value in sorted(counters.items(), key=lambda x: (
        not x[0].startswith("TOTAL"), "genotype" in x[0], "variants" not in x[0], len(x[0]), -1*x[1])
):
    if "TOTAL" in key:
        continue

    denominator_key = "TOTAL variants" if "genotype" in key or "variants" in key else "TOTAL alleles"
    total = counters[denominator_key]
    percent = 100 * value / total

    if percent < args.min_percent:
        continue
    
    try:
        percent = f"{percent:5.1f}%"
    except:
        percent = ""

    print(f"{value:10,d} ({percent}) {key:17s}  (out of {total:9,d} {denominator_key})")

