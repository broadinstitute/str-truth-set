"""This script takes a VCF file and prints summary stats about the number of variants and alleles"""

import argparse
import collections
import gzip

p = argparse.ArgumentParser()
p.add_argument("-p", "--prefix", help="Optional prefix string to print at the beginning of every line")
p.add_argument("-m", "--min-percent", default=0, type=float,
               help="Only print stats (eg. % of variants that is multi-allelic) if the % is greater than this threshold")
p.add_argument("-l", "--label", help="Optional text label to describe the VCF")
p.add_argument("--count-by-filter", action="store_true", help="Whether to print stats separately by filter field")
p.add_argument("vcf_path")
args = p.parse_args()

prefix = ""
if args.prefix:
    prefix = f"{args.prefix}      "

if args.label:
    print(f"{prefix}VCF stats for {args.label}",  "-- "*50, args.vcf_path)
else:
    print(f"{prefix}VCF stats for ",  "-- "*50, args.vcf_path)

# step 1: Compute stats
counters = collections.defaultdict(int)
fopen = gzip.open if args.vcf_path.endswith("gz") else open
for i, line in enumerate(fopen(args.vcf_path, "rt")):
    if line.startswith("#"):
        continue

    # VCF Header:  #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  syndip
    fields = line.strip().split("\t")
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
        if alt == "*":
            continue

        counters["TOTAL alleles"] += 1

        if len(ref) == len(alt):
            snv_counter += 1
            if len(ref) > 1:
                counters[f"MNV alleles"] += 1
            else:
                counters[f"SNV alleles"] += 1
        elif len(ref) > len(alt) and ref.startswith(alt):
            deletion_counter += 1
            counters[f"DEL alleles"] += 1
        elif len(ref) < len(alt) and alt.startswith(ref) == 1:
            insertion_counter += 1
            counters[f"INS alleles"] += 1
        else:
            counters[f"complex MNV indel alleles"] += 1

    # count variant types
    variant_filter_prefixes = ["", f"filter:{fltr}  ", "filtered:  "] if args.count_by_filter else [""]
    for variant_filter_prefix in variant_filter_prefixes:
        if insertion_counter > 0 and deletion_counter == 0 and snv_counter == 0:
            if insertion_counter > 1:
                counters[f"{variant_filter_prefix}multiallelic INS variants"] += 1
            else:
                counters[f"{variant_filter_prefix}INS variants"] += 1
        elif insertion_counter == 0 and deletion_counter > 0 and snv_counter == 0:
            if deletion_counter > 1:
                counters[f"{variant_filter_prefix}multiallelic DEL variants"] += 1
            else:
                counters[f"{variant_filter_prefix}DEL variants"] += 1
        elif insertion_counter == 0 and deletion_counter == 0 and snv_counter > 0:
            if snv_counter > 1:
                counters[f"{variant_filter_prefix}multiallelic SNV variants"] += 1
            else:
                counters[f"{variant_filter_prefix}SNV variants"] += 1
        else:
            var_types = []
            if insertion_counter > 0:
                var_types.append("INS")
            if deletion_counter > 0:
                var_types.append("DEL")
            if snv_counter > 0:
                var_types.append("SNV")

            var_types = "/".join(var_types)
            counters[f"{variant_filter_prefix}mixed multiallelic {var_types} variants"] += 1

# step 2: Print totals
print(prefix+"{:10,d}".format(counters["TOTAL variants"]), "TOTAL variants")
print(prefix+"{:10,d}".format(counters["TOTAL alleles"]), "TOTAL alleles",
    ("(%0.2f alleles per variant)" % (counters["TOTAL alleles"]/counters["TOTAL variants"])) if counters["TOTAL variants"] > 0 else ""
)

# step 3: Print all counters
for keyword in "genotype", "variants", "alleles":
    current_counter = [(key, count) for key, count in counters.items() if keyword in key]
    if args.count_by_filter:
        sort_key = lambda x: (x[0], -x[1])
    else:
        sort_key = lambda x: (-x[1], x[0])


    current_counter = sorted(current_counter, key=sort_key)
    print("--------------")
    for key, value in current_counter:
        if "variants" in key or "genotype" in key:
            total_key = "TOTAL variants"
        else:
            total_key = "TOTAL alleles"

        total = counters[total_key]
        percent = f"{100*value / total:5.1f}%" if total > 0 else ""

        print(f"{prefix}{value:10,d} out of {total:10,d} ({percent}) {key}")

print("--------------")
