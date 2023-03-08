import collections
import pandas as pd
from figures_and_tables.numbers_utils import format_np

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat == "Yes"]
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_alleles = df_alleles[df_alleles.IsPureRepeat == "Yes"]

#print("\n".join(df_variants.columns))
columns = ["OverlapsIlluminaSTRCatalog: Locus",
           "OverlapsGangSTRCatalog13: Locus", "OverlapsGangSTRCatalog17: Locus",
           "OverlapsHipSTRCatalog: Locus"] + [f"OverlapsTRFPureRepeats{i}bp: Locus" for i in (6, 9, 12, 15)]

#%%

print("    -----")
print("% STR truth set loci missed by different catalogs: ")
for column in columns:
    count_values = collections.Counter(df_variants[column])
    print(format_np(count_values["no overlap"], len(df_variants)), "of STR variants don't overlap", column.replace("Overlaps", "").split(":")[0])

#%%

print("    -----")
print("% STR truth set loci missed by different catalogs for loci with a reference locus >= 24bp: ")

df_variants_filtered = df_variants[df_variants.NumRepeatsInReference * df_variants.MotifSize >= 24]
print(format_np(len(df_variants_filtered), len(df_variants)), "have a reference locus size of at least 24bp")
for column in columns:
    count_values = collections.Counter(df_variants_filtered[column])
    print(format_np(count_values["no overlap"], len(df_variants_filtered)), "of STR variants don't overlap", column.replace("Overlaps", "").split(":")[0])


#%%

print("    -----")
print("% STR truth set loci missed by different catalogs for loci with a reference locus >= 5 repeats: ")

df_variants_filtered = df_variants[(df_variants.NumRepeatsLongAllele - df_variants.NumRepeatsInReference) >= 5]
print(format_np(len(df_variants_filtered), len(df_variants)), "have a reference locus >= 5 repeats")
for column in columns:
    count_values = collections.Counter(df_variants_filtered[column])
    print(format_np(count_values["no overlap"], len(df_variants_filtered)), "of STR variants don't overlap", column.replace("Overlaps", "").split(":")[0])

#%%

print("    -----")
print("% STR truth set loci missed by different catalogs for loci with a reference locus >= 10 repeats: ")

df_variants_filtered = df_variants[(df_variants.NumRepeatsLongAllele - df_variants.NumRepeatsInReference) >= 10]
print(format_np(len(df_variants_filtered), len(df_variants)), "have a reference locus >= 10 repeats")
for column in columns:
    count_values = collections.Counter(df_variants_filtered[column])
    print(format_np(count_values["no overlap"], len(df_variants_filtered)), "of STR variants don't overlap", column.replace("Overlaps", "").split(":")[0])
