import collections
import gzip
import pandas as pd

print_variant_counts = False

df_TR_alleles = pd.read_table("step2.STRs.alleles.tsv.gz")

f = gzip.open("step2.STRs.filtered_out_indels.vcf.gz", "rt")
variant_counters = collections.defaultdict(int)
allele_counters = collections.defaultdict(int)

# numbers from paper_numbers1_results_part1_defining_the_truth_set.py
total_indel_variants = 518_285
total_indel_alleles = 535_354

allele_counters["total indel alleles"] = total_indel_alleles
allele_counters["total TR alleles in truth set"] = len(df_TR_alleles)

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

allele_counters["sum"] = sum(allele_counters[k] for k in allele_counters.keys() if k != "total indel alleles")

print("Allele counts:")
for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
    print(f"{count:10,d} ({100 * count / total_indel_alleles:5.1f}%) {key:10s}")

print("HTML table:")
print("<table>")
print(f"<tr><th>Description</th><th># of indels</th><th>% of indels</th></tr>")
for i, (key, count) in enumerate(sorted(allele_counters.items(), key=lambda x: -x[1])):
    print(f"<tr><td>{key}</td><td>{count:10,d}</td><td>{100 * count / total_indel_alleles:5.1f}%</td></tr>")
    if i == 0:
        print(f"<tr><td></td><td></td><td></td></tr>")
print("</table>")

print("------")

if print_variant_counts:
    for key, count in sorted(variant_counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} ({100 * count / total_indel_variants:5.1f}%) {key:10s}")

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

        print(f"{value:10,d} ({100 * value / total_indel_variants:5.1f}%) {label:10s}")



#%%