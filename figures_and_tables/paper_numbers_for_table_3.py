import collections
import gzip
import pandas as pd

print_variant_counts = False

df = pd.read_table("step2.STRs.alleles.tsv.gz")

total_STR_alleles_in_STR_vcf = 0
f1 = gzip.open("step2.STRs.vcf.gz", "rt")
f1 = gzip.open("temp/step2.pure_STRs.vcf.gz", "rt")
for line in f1:
    if line.startswith("#"):
        continue
    fields = line.split("\t")
    alt_alleles = fields[4].split(",")
    total_STR_alleles_in_STR_vcf += len(alt_alleles)

f = gzip.open("temp/step2.pure_STRs.filtered_out_indels.vcf.gz", "rt")
variant_counters = collections.defaultdict(int)
allele_counters = collections.defaultdict(int)

# numbers from paper_numbers1_results_part1_defining_the_truth_set.py
total_variants = 518_285
total_alleles = 286_446 + 271_635
total_alleles = 285_701 + 271_124

allele_counters["total alleles"] = total_alleles
allele_counters["total TRs in truth set"] = total_STR_alleles_in_STR_vcf #len(df)

for line in f:
    if line.startswith("#"):
        continue
    fields = line.split("\t")
    filter_values = fields[6]
    variant_counters[filter_values] += 1
    for filter_value in filter_values.split(";"):
        if filter_value in ("SNV/MNV", "complex multinucleotide insertion + deletion"):
            continue
        allele_counters[filter_value] += 1

for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
    print(f"<tr><td>{count:10,d} ({100 * count / total_variants:5.1f}%) {key:10s}")

allele_counters["sum"] = sum(allele_counters[k] for k in allele_counters.keys() if k != "total alleles")

for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
    print(f"{count:10,d} ({100 * count / total_alleles:5.1f}%) {key:10s}")

print("total_alleles_in_STR_vcf; ", total_STR_alleles_in_STR_vcf)
print("------")

if print_variant_counts:
    for key, count in sorted(variant_counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} ({100 * count / total_variants:5.1f}%) {key:10s}")

    print("---")
    other_counter_keys = set(variant_counters.keys())
    for label, counter_keys in {
        "homopolymer": ("repeat unit < 2 bp", "repeat unit < 2 bp;repeat unit < 2 bp",),
        "spans < 9 bp": ("spans < 9 bp", "spans < 9 bp;spans < 9 bp",),
        "INDEL without repeats": ("INDEL without repeats", "INDEL without repeats;INDEL without repeats",),
        "INDEL with only 2 repeats": ("is only 2 repeats", "is only 2 repeats;is only 2 repeats",),
        "motif > 50bp": ("repeat unit > 50 bp", "repeat unit > 50 bp;repeat unit > 50 bp",),
        "ends in partial repeats": ("ends in partial repeat",),
        "multi-allelic with SNV": ("SNV/MNV;STR allele", "SNV/MNV;repeat unit < 2 bp", "STR allele;SNV/MNV", "SNV/MNV;repeat unit > 50 bp",),
        #"locus overlaps more than one STR variant": ("locus overlaps more than one STR variant",),
        "other": ("other",),
    }.items():
        if counter_keys == ("other",):
            value = sum(variant_counters[key] for key in other_counter_keys)
        else:
            value = sum(variant_counters[key] for key in counter_keys)
            other_counter_keys = other_counter_keys - set(counter_keys)

        print(f"{value:10,d} ({100 * value / total_variants:5.1f}%) {label:10s}")



#%%