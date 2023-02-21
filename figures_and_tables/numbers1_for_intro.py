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

#%%

# Total STR stats

print("========")
print("Truth set: pure STR loci\n")

print(f"{format_n(len(df_variants))}    TOTAL STR loci")
print(f"{format_n(len(df_alleles))}    TOTAL STR alleles")

#%%

with open("step1.log", "rt") as f:
    step1_log_contents = f.read()

#%%
print("-"*100)
print("Numbers for intro:")
print("-"*100)
percent_missed_by_gangstr_catalog = re.search(
    f"GangSTRCatalog17.*{len(df_variants):,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    step1_log_contents).group(1)
print(f"{percent_missed_by_gangstr_catalog}% missed by GangSTR catalog")

percent_missed_by_illumina_catalog = re.search(
    f"IlluminaSTRCatalog.*{len(df_variants):,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    step1_log_contents).group(1)
print(f"{percent_missed_by_illumina_catalog}% missed by ExpansionHunter catalog")

#%%

print("-"*100)
# STR alleles > 30bp relative to reference
bp_diff_from_ref = (df_alleles["RepeatSize (bp)"] - df_alleles["NumRepeatsInReference"] * df_alleles["MotifSize"]).astype(int).abs()

assert all(bp_diff_from_ref > 1)

print(f"{format_np(sum(bp_diff_from_ref > 30), len(df_alleles))}    STR alleles that differ from reference locus size by more than 30bp")


#%%

print("-"*100)

# Novel STR stats
is_found_in_reference_counter = collections.Counter(df_variants.IsFoundInReference)
print(f"{format_np(is_found_in_reference_counter['Yes'], len(df_variants))}    FOUND IN REF STR loci (found in reference genome)")
print(f"{format_np(is_found_in_reference_counter['No'], len(df_variants))}    NOVEL STR loci (not found in reference genome)")

assert is_found_in_reference_counter['Yes'] + is_found_in_reference_counter['No'] == len(df_variants)
#%%

print("-"*100)
print("Defining the STR truth set:")
print("-"*100)
# Defining the STR truth set - numbers for figure
total_variants = int(re.search(
    f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))
total_high_confidence_variants = int(re.search(
    f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

total_alleles = int(re.search(
    f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))
total_high_confidence_alleles = int(re.search(
    f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))

high_confidence_INS_alleles = int(re.search(
    f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*INS alleles", step1_log_contents).group(1).replace(",", ""))
high_confidence_DEL_alleles = int(re.search(
    f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*DEL alleles", step1_log_contents).group(1).replace(",", ""))


high_confidence_SNV_variants = int(re.search(
    f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_variants:,d}.*SNV variants", step1_log_contents).group(1).replace(",", ""))

high_confidence_INS_DEL_variants = total_high_confidence_variants - high_confidence_SNV_variants

total_STR_variants_before_validation_step = int(re.search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))
total_STR_alleles_before_validation_step = int(re.search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))


print(f"{format_n(total_variants)} total variants in raw SynDip")
print(f"{format_n(total_alleles)} total alleles in raw SynDip")

print(f"{format_n(total_high_confidence_variants)} high-confidence SynDip variants")
print(f"{format_n(total_high_confidence_alleles)} high-confidence SynDip alleles")

print(f"{format_n(high_confidence_INS_DEL_variants)} high-confidence INS or DEL variants")

print(f"{format_np(high_confidence_INS_alleles, total_high_confidence_alleles)} high-confidence INS alleles")
print(f"{format_np(high_confidence_DEL_alleles, total_high_confidence_alleles)} high-confidence DEL alleles")
print(f"{format_np(high_confidence_INS_alleles + high_confidence_DEL_alleles, total_high_confidence_alleles)} high-confidence INS + DEL alleles")

print(f"{format_n(total_STR_variants_before_validation_step)} total pure STR variants")
print(f"{format_np(total_STR_alleles_before_validation_step, total_high_confidence_alleles)} pure STR alleles")

