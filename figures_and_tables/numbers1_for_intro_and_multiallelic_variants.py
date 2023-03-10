"""This script computes numbers that appear in the text of the paper, as well as all tables"""

import collections
import pandas as pd
from figures_and_tables.numbers_utils import format_n, format_np, search

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat == "Yes"]
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_alleles = df_alleles[df_alleles.IsPureRepeat == "Yes"]

#%%

# Total STR stats
print("========")
print("Truth set: pure STR loci\n")

print(f"{format_n(len(df_variants))}    TOTAL STR loci")
print(f"{format_n(len(df_alleles))}    TOTAL STR alleles")

#%%

with open("step_A.log", "rt") as f:
    stepA_log_contents = f.read()

#%%

print("-"*100)
# STR alleles > 30bp relative to reference
bp_diff_from_ref = (df_alleles["RepeatSize (bp)"] - df_alleles["NumRepeatsInReference"] * df_alleles["MotifSize"]).astype(int).abs()

assert all(bp_diff_from_ref > 1)

print(f"{format_np(sum(bp_diff_from_ref > 30), len(df_alleles))}    STR alleles that differ from reference locus size by more than 30bp")

print(f"{format_np(sum(bp_diff_from_ref > 60), len(df_alleles))}    STR alleles that differ from reference locus size by more than 60bp")

#%%

print("-"*100)
print("Numbers for intro:")
print("-"*100)
total = len(df_variants[df_variants.IsFoundInReference == "Yes"])
percent_missed_by_gangstr_catalog = search(
    f"GangSTRCatalog17.*{total:,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    stepA_log_contents)
print(f"{percent_missed_by_gangstr_catalog}% missed by GangSTR catalog")

percent_missed_by_illumina_catalog = search(
    f"IlluminaSTRCatalog.*{total:,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    stepA_log_contents)
print(f"{percent_missed_by_illumina_catalog}% missed by ExpansionHunter catalog")

percent_missed_by_hipstr_catalog = search(
    f"HipSTRCatalog.*{total:,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    stepA_log_contents)
print(f"{percent_missed_by_hipstr_catalog}% missed by HipSTR catalog")

percent_missed_by_9bp_catalog = search(
    f"TRFPureRepeats9bp.*{total:,d} [(][ ]+([0-9]+[.][0-9]+)[%][)].* truth set loci:.*no overlap",
    stepA_log_contents)
print(f"{percent_missed_by_9bp_catalog}% missed by TRFPureRepeats9bp catalog")

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
total_variants = int(search(f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))

total_high_confidence_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))
total_high_confidence_monoallelic_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] SNV variants", stepA_log_contents).replace(",", ""))
total_high_confidence_monoallelic_INS_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] INS variants", stepA_log_contents).replace(",", ""))
total_high_confidence_monoallelic_DEL_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] DEL variants", stepA_log_contents).replace(",", ""))

total_high_confidence_multiallelic_variants = total_high_confidence_variants - (
    total_high_confidence_monoallelic_SNV_variants +
    total_high_confidence_monoallelic_INS_variants +
    total_high_confidence_monoallelic_DEL_variants)


total_alleles = int(search(f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents).replace(",", ""))
total_high_confidence_alleles = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents).replace(",", ""))

high_confidence_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_variants:,d}.*[)] SNV variants", stepA_log_contents).replace(",", ""))
high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))

total_STR_variants_before_validation_step = int(search(f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))
total_STR_alleles_before_validation_step = int(search(f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents).replace(",", ""))


total_high_confidence_multiallelic_INS_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] multiallelic INS variants", stepA_log_contents).replace(",", ""))
total_high_confidence_multiallelic_DEL_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] multiallelic DEL variants", stepA_log_contents).replace(",", ""))
total_high_confidence_multiallelic_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] multiallelic SNV variants", stepA_log_contents).replace(",", ""))
total_high_confidence_mixed_multiallelic_INS_DEL_variants = int(search(f"step1:output:[ ]*([0-9,]+)[ ]* out of[ ]* {total_high_confidence_variants:,d}.*[)] mixed multiallelic INS/DEL variants", stepA_log_contents).replace(",", ""))


total_high_confidence_SNV_variants = (
        total_high_confidence_monoallelic_SNV_variants +
        total_high_confidence_multiallelic_SNV_variants)

total_high_confidence_INS_variants = (
        total_high_confidence_monoallelic_INS_variants +
        total_high_confidence_multiallelic_INS_variants +
        total_high_confidence_mixed_multiallelic_INS_DEL_variants)
total_high_confidence_DEL_variants = (
        total_high_confidence_monoallelic_DEL_variants +
        total_high_confidence_multiallelic_DEL_variants +
        total_high_confidence_mixed_multiallelic_INS_DEL_variants)

total_high_confidence_multiallelic_INDEL_variants = (
        total_high_confidence_multiallelic_INS_variants +
        total_high_confidence_multiallelic_DEL_variants +
        total_high_confidence_mixed_multiallelic_INS_DEL_variants)


total_high_confidence_INDEL_variants = total_high_confidence_variants - total_high_confidence_SNV_variants

print(f"{format_n(total_variants)} total variants in raw SynDip")
print(f"{format_n(total_alleles)} total alleles in raw SynDip")

print(f"{format_n(total_high_confidence_variants)} high-confidence SynDip variants")
print(f"{format_n(total_high_confidence_alleles)} high-confidence SynDip alleles")

print(f"{format_n(total_high_confidence_INDEL_variants)} high-confidence INS or DEL variants")
print(f"{format_np(total_high_confidence_multiallelic_variants, total_high_confidence_variants)} high confidence multiallelic variants")
print(f"{format_np(high_confidence_INS_alleles, total_high_confidence_alleles)} high-confidence INS alleles")
print(f"{format_np(high_confidence_DEL_alleles, total_high_confidence_alleles)} high-confidence DEL alleles")
print(f"{format_np(high_confidence_INS_alleles + high_confidence_DEL_alleles, total_high_confidence_alleles)} high-confidence INS + DEL alleles")

print(f"{format_n(total_STR_variants_before_validation_step)} total pure STR variants")
print(f"{format_np(total_STR_alleles_before_validation_step, total_high_confidence_alleles)} pure STR alleles")


#%%

is_multiallelic_counter = collections.Counter(df_variants.IsMultiallelic)
total_multiallelic_STR_variants = is_multiallelic_counter["Yes"]

print("Fraction of STRs that is multi-allelic: ")
print(f"{format_np(total_multiallelic_STR_variants, len(df_variants))} fraction multiallelic STR variants")

print("Fraction of other variants that is multi-allelic: ")
print(f"{format_np(total_high_confidence_multiallelic_SNV_variants, total_high_confidence_SNV_variants)} fraction high-confidence multi-allelic SNVs in SynDip")
print(f"{format_np(total_high_confidence_multiallelic_INS_variants + total_high_confidence_mixed_multiallelic_INS_DEL_variants, total_high_confidence_INS_variants)} fraction high-confidence multi-allelic INS in SynDip")
print(f"{format_np(total_high_confidence_multiallelic_DEL_variants + total_high_confidence_mixed_multiallelic_INS_DEL_variants, total_high_confidence_DEL_variants)} fraction high-confidence multi-allelic DEL in SynDip")
print(f"{format_np(total_high_confidence_multiallelic_INDEL_variants, total_high_confidence_INDEL_variants)} fraction high-confidence multi-allelic INS or DEL in SynDip")
print(f"{format_np(total_high_confidence_multiallelic_INDEL_variants - total_multiallelic_STR_variants, total_high_confidence_INDEL_variants - len(df_variants))} fraction high-confidence multi-allelic INS or DEL in SynDip minus STRs")

#%%


