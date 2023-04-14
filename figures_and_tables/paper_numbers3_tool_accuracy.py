import pandas as pd
from figures_and_tables.numbers_utils import format_np

#%%
df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
#df_variants = df_variants[df_variants.IsPureRepeat]
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
#df_alleles = df_alleles[df_alleles.IsPureRepeat]

df = pd.read_table("./tool_comparison/combined.results.variants.tsv.gz")
df = df[(df.Coverage == "30x") & (df.PositiveOrNegative == "positive")]


#%%

df_2to6bp = df[(df["MotifSize"] >= 2) & (df["MotifSize"] <= 6)]

print(f"{format_np(len(df_2to6bp[df_2to6bp['NumRepeats: Allele 2: HipSTR'].isna()]), len(df_2to6bp))} of variants were not called by HipSTR")

print(f"{format_np(len(df_2to6bp[~df_2to6bp['NumRepeats: Allele 2: ExpansionHunter'].isna() & ~df_2to6bp['NumRepeats: Allele 2: GangSTR'].isna() & ~df_2to6bp['NumRepeats: Allele 2: HipSTR'].isna()]), len(df_2to6bp))}  of genotypes were called by all 3 tools")


#%%

df_variants_vs_EHdn = pd.read_table(
    "./tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/"
    "CHM1_CHM13_WGS2.downsampled_to_30x.truth_set_EHdn_comparison_table.tsv")

print(list(df_variants_vs_EHdn.columns))

import collections
collections.Counter(df_variants_vs_EHdn["EHdn Concordance"])

df_variants_over_150bp = df_variants_vs_EHdn[(df_variants_vs_EHdn["RepeatSizeLongAllele (bp)"] >= 150)]
df_pure_variants_over_150bp = df_variants_over_150bp[df_variants_over_150bp.IsPureRepeat]

df_2to6bp = df_pure_variants_over_150bp[df_pure_variants_over_150bp["MotifSize"] <= 6]
df_7to24bp = df_pure_variants_over_150bp[(df_pure_variants_over_150bp["MotifSize"] <= 24) & (df_pure_variants_over_150bp["MotifSize"] >= 7)]
df_25to50bp = df_pure_variants_over_150bp[(df_pure_variants_over_150bp["MotifSize"] <= 50) & (df_pure_variants_over_150bp["MotifSize"] >= 25)]

print(f"{format_np(sum(df_2to6bp['EHdn Concordance'] != 'No Call'), len(df_2to6bp))} of all pure 2-6bp variants over 150bp were called by EHdn")
print(f"{format_np(sum(df_7to24bp['EHdn Concordance'] != 'No Call'), len(df_7to24bp))} of all pure 7-24bp variants over 150bp were called by EHdn")
print(f"{format_np(sum(df_25to50bp['EHdn Concordance'] != 'No Call'), len(df_25to50bp))} of all pure 25-50bp variants over 150bp were called by EHdn")

print(f"{format_np(sum(df_pure_variants_over_150bp['EHdn Concordance'] != 'No Call'), len(df_pure_variants_over_150bp))} of all pure variants over 150bp were called by EHdn")

#%%

df_EHdn_vs_variants = pd.read_table("./tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/CHM1_CHM13_WGS2.downsampled_to_30x.EHdn_results_table.with_truth_set_concordance.tsv")
df_EHdn_vs_variants_2to6bp = df_EHdn_vs_variants[(df_EHdn_vs_variants["MotifSize"] >= 2) & (df_EHdn_vs_variants["MotifSize"] <= 6)]
df_EHdn_vs_variants_7to24bp = df_EHdn_vs_variants[(df_EHdn_vs_variants["MotifSize"] >= 7) & (df_EHdn_vs_variants["MotifSize"] <= 24)]
df_EHdn_vs_variants_25to50bp = df_EHdn_vs_variants[(df_EHdn_vs_variants["MotifSize"] >= 25) & (df_EHdn_vs_variants["MotifSize"] <= 50)]

print(f"{format_np(sum(df_EHdn_vs_variants['EHdn Concordance With Truth Set'] == 'False Positive'), len(df_EHdn_vs_variants))} false positives")

print(f"{format_np(sum(df_EHdn_vs_variants_2to6bp['EHdn Concordance With Truth Set'] == 'False Positive'), len(df_EHdn_vs_variants_2to6bp))} false positives for 2-6bp variants")
print(f"{format_np(sum(df_EHdn_vs_variants_7to24bp['EHdn Concordance With Truth Set'] == 'False Positive'), len(df_EHdn_vs_variants_7to24bp))} false positives for 7-24bp variants")
print(f"{format_np(sum(df_EHdn_vs_variants_25to50bp['EHdn Concordance With Truth Set'] == 'False Positive'), len(df_EHdn_vs_variants_25to50bp))} false positives for 25-50bp variants")

#%%


#%%
