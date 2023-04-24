import pandas as pd
from figures_and_tables.numbers_utils import format_np

#%%
df_before_validation = pd.read_table("step2.STRs.variants.tsv.gz")
df_raw_before_validation = pd.read_table("step2.raw_STRs_including_homopolymers.variants.tsv.gz")
df_raw_before_validation = df_raw_before_validation[(df_raw_before_validation.MotifSize > 1) & (df_raw_before_validation.MotifSize <= 50)]

df_after_validation = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
#df_raw_after_validation = pd.read_table("raw_STRs_including_homopolymer_truth_set.v1.variants.tsv.gz")
df_raw_after_validation = pd.read_table("step7.filtered_raw_STRs_including_homopolymers.variants.tsv.gz")
df_raw_after_validation = df_raw_after_validation[(df_raw_after_validation.MotifSize > 1) & (df_raw_after_validation.MotifSize <= 50)]

print(f"{format_np(len(df_raw_after_validation), len(df_raw_before_validation))} raw TRs would have passed all validation steps.")
print(f"Which is {format_np(len(df_raw_after_validation) - len(df_after_validation), len(df_raw_before_validation) - len(df_before_validation))} additional TRs relative to the truth set.")

#%%