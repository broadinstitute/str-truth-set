"""This script computes numbers that appear in the text of the paper, as well as all tables"""

import collections
import pandas as pd
import re

df_variants = pd.read_table("STR_truthset.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat == "Yes"]
df_alleles = pd.read_table("STR_truthset.v1.alleles.tsv.gz")
df_alleles = df_alleles[df_alleles.IsPureRepeat == "Yes"]

#%%


def format_p(count, total):
    return f"{100*count/total:5.1f}%"


def format_n(count, d=10):
    return f"{count:{d},d}"


def format_np(count, total, d=10):
    return f"{format_n(count, d=d)} out of {total:{d},d}  ({format_p(count, total)})"

# Totals

print("========")
print("Truth set: pure STR loci\n")

print(f"{format_n(len(df_variants))}    TOTAL STR loci")
print(f"{format_n(len(df_alleles))}    TOTAL STR alleles")

#%%

#146,618 ( 41.3%) truth set loci:     no overlap

with open("step1.log", "rt") as f:
    step1_log_contents = f.read()
#%%
percent_missed_by_gangstr_catalog = re.search(
    f"GangSTRCatalog17.*{len(df_variants):,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    step1_log_contents).group(1)
print(f"{percent_missed_by_gangstr_catalog}% missed by GangSTR catalog")

percent_missed_by_illumina_catalog = re.search(
    f"IlluminaSTRCatalog.*{len(df_variants):,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    step1_log_contents).group(1)
print(f"{percent_missed_by_illumina_catalog}% missed by ExpansionHunter catalog")

#%%

# STR alleles > 30bp relative to reference

bp_diff_from_ref = (df_alleles["RepeatSize (bp)"] - df_alleles["NumRepeatsInReference"] * df_alleles["MotifSize"]).astype(int).abs()

assert all(bp_diff_from_ref > 1)

print(f"{format_np(sum(bp_diff_from_ref > 30), len(df_alleles))}    STR alleles that differ from reference locus size by more than 30bp")


#%%

# Novel STRs

is_found_in_reference_counter = collections.Counter(df_variants.IsFoundInReference)
print(f"{format_np(is_found_in_reference_counter['Yes'], len(df_variants))}    FOUND IN REF STR loci (found in reference genome)")
print(f"{format_np(is_found_in_reference_counter['No'], len(df_variants))}    NOVEL STR loci (not found in reference genome)")

assert is_found_in_reference_counter['Yes'] + is_found_in_reference_counter['No'] == len(df_variants)
#%%


#%%