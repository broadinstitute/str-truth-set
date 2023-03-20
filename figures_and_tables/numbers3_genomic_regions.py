import collections
import pandas as pd
from figures_and_tables.numbers_utils import format_np

df_variants = pd.read_table("pure_STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat]
df_alleles = pd.read_table("pure_STR_truth_set.v1.alleles.tsv.gz")
df_alleles = df_alleles[df_alleles.IsPureRepeat]

#%%

num_STRs_overlap_segdups = sum(df_variants.OverlapsSegDupIntervals)
print(f"{format_np(num_STRs_overlap_segdups, len(df_variants))}    STR loci overlap segdups")

for gene_region_column in "GeneRegionFromGencode_V42", "GeneRegionFromMane_V1":
    print("-"*100)
    print(gene_region_column)
    for gene_region, count in sorted(collections.Counter(df_variants[gene_region_column]).items(), key=lambda i:i[0]):
        print(format_np(count, len(df_variants)), f" {gene_region} from {gene_region_column}")

#%%


df_MANE_CDS_alleles = df_alleles[df_alleles["GeneRegionFromMane_V1"] == "CDS"]

print(format_np(sum(df_MANE_CDS_alleles.MotifSize % 3 == 0), len(df_MANE_CDS_alleles)), " of coding STR alleles in MANE v1 have a motif size that's a multiple of 3bp")
print("Coding STR alleles where the motif size isn't a multiple of 3bp (and so could be an out-of-frame insertion/deletion)")
df_MANE_CDS_alleles_not_multiple_of_3bp = df_MANE_CDS_alleles[df_MANE_CDS_alleles.MotifSize % 3 != 0]
print("\n".join(df_MANE_CDS_alleles_not_multiple_of_3bp.SummaryString))
print(df_MANE_CDS_alleles_not_multiple_of_3bp[["Locus", "MotifSize", "Motif", "RepeatSize (bp)", "NumRepeats"]])

#%%

print(format_np(len(df_variants[~df_variants.IsFoundInReference]), len(df_variants)), "are novel STR loci")
print(format_np(len(df_variants[df_variants.IsFoundInReference]), len(df_variants)), "of STRs exist in the hg38")
print("Gene regions of novel variants: ", set(df_variants[~df_variants.IsFoundInReference].GeneRegionFromGencode_V42))
#%%


#%%