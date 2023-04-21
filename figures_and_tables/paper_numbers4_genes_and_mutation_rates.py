import pandas as pd
from figures_and_tables.numbers_utils import format_np

#%%
df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")


#%%

non_exonic_variants = df_variants[df_variants.GeneRegionFromMane_V1.isin({'intergenic', 'intron'})]
print(f"{format_np(len(non_exonic_variants), len(df_variants))} "
      f"variants in the truth set are within MANE v1 intergenic or intronic regions")

coding_variants_with_pure_repeats = df_variants[(df_variants.GeneRegionFromMane_V1 == "CDS") & df_variants.IsPureRepeat]

print(f"{format_np(sum(coding_variants_with_pure_repeats.MotifSize % 3 == 0), len(coding_variants_with_pure_repeats))} "
      f"variants in the truth set overlap MANE v1 coding regions")

#%%
#num_MANE_genes = len(set(coding_variants_with_pure_repeats.GeneNameFromMane_V1))
#print(f"{format_np(len(df_alleles), num_MANE_genes)} "
#coding_variants_with_pure_repeats