import pandas as pd
from figures_and_tables.numbers_utils import format_np

#%%
df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
#df_variants = df_variants[df_variants.IsPureRepeat]
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
#df_alleles = df_alleles[df_alleles.IsPureRepeat]

df = pd.read_table("./tool_comparison/combined.results.alleles.tsv.gz")
df = df[(df.Coverage == "40x") & (df.PositiveOrNegative == "positive")]


#%%

print(f"{format_np(len(df[df['NumRepeats: Allele: HipSTR'].isna()]), len(df))} of alleles were not called by HipSTR")

print(f"{format_np(len(df[~df['NumRepeats: Allele: ExpansionHunter'].isna() & ~df['NumRepeats: Allele: GangSTR'].isna() & ~df['NumRepeats: Allele: HipSTR'].isna()]), len(df))}  of genotypes were called by all 3 tools")


#%%