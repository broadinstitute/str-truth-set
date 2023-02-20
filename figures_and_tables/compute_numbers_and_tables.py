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

# Novel STR stats

is_found_in_reference_counter = collections.Counter(df_variants.IsFoundInReference)
print(f"{format_np(is_found_in_reference_counter['Yes'], len(df_variants))}    FOUND IN REF STR loci (found in reference genome)")
print(f"{format_np(is_found_in_reference_counter['No'], len(df_variants))}    NOVEL STR loci (not found in reference genome)")

assert is_found_in_reference_counter['Yes'] + is_found_in_reference_counter['No'] == len(df_variants)
#%%

# Defining the STR truth set - numbers for figure
total_variants = int(re.search(
    f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))
total_high_confidence_variants = int(re.search(
    f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

total_alleles = int(re.search(
    f"step1:input:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))
total_high_confidence_alleles = int(re.search(
    f"step1:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))

high_confidence_SNV_variants = int(re.search(
    f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_variants:,d}.*SNV variants", step1_log_contents).group(1).replace(",", ""))

high_confidence_INS_DEL_variants = total_high_confidence_variants - high_confidence_SNV_variants
high_confidence_INS_DEL_alleles = (
    int(re.search(
        f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*INS alleles", step1_log_contents).group(1).replace(",", ""))
    + int(re.search(
        f"step1:output:[ ]*([0-9,]+) out of[ ]* {total_high_confidence_alleles:,d}.*DEL alleles", step1_log_contents).group(1).replace(",", "")))


total_STR_variants_before_validation_step = int(re.search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))
total_STR_alleles_before_validation_step = int(re.search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", step1_log_contents).group(1).replace(",", ""))


print(f"{format_n(total_variants)} total variants in raw SynDip")
print(f"{format_n(total_alleles)} total alleles in raw SynDip")

print(f"{format_n(total_high_confidence_variants)} high-confidence SynDip variants")
print(f"{format_n(total_high_confidence_alleles)} high-confidence SynDip alleles")

print(f"{format_n(high_confidence_INS_DEL_variants)} high-confidence INS or DEL variants")
print(f"{format_np(high_confidence_INS_DEL_alleles, total_high_confidence_alleles)} high-confidence INS or DEL alleles")

print(f"{format_n(total_STR_variants_before_validation_step)} total pure STR variants")
print(f"{format_np(total_STR_alleles_before_validation_step, total_high_confidence_alleles)} pure STR alleles")

#%%

# STR Validation Stats - 1st row of figure:

print(f"{format_n(total_STR_variants_before_validation_step)}    TOTAL STR loci before validation")

def look_up_liftover_stats(prefix, total):
    INS_STRs_before_validation_step = int(re.search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.* INS variants", step1_log_contents).group(1).replace(",", ""))

    DEL_STRs_before_validation_step = int(re.search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.* DEL variants", step1_log_contents).group(1).replace(",", ""))

    multiallelic_INS_STRs_before_validation_step = int(re.search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.* multiallelic INS", step1_log_contents).group(1).replace(",", ""))
    multiallelic_DEL_STRs_before_validation_step = int(re.search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.* multiallelic DEL", step1_log_contents).group(1).replace(",", ""))
    mixed_multiallelic_STRs_before_validation_step = int(re.search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.* mixed multiallelic", step1_log_contents).group(1).replace(",", ""))

    assert total == (
            INS_STRs_before_validation_step + DEL_STRs_before_validation_step + multiallelic_DEL_STRs_before_validation_step
            + multiallelic_INS_STRs_before_validation_step + mixed_multiallelic_STRs_before_validation_step)

    print(f"{format_np(INS_STRs_before_validation_step + multiallelic_INS_STRs_before_validation_step, total)} "
          f"INS STRs in {prefix}")
    print(f"{format_np(multiallelic_INS_STRs_before_validation_step, INS_STRs_before_validation_step + multiallelic_INS_STRs_before_validation_step)} "
          f"are multi-allelic INS")
    print(f"{format_np(DEL_STRs_before_validation_step + multiallelic_DEL_STRs_before_validation_step, total)} "
          f"DEL STRs in {prefix}")
    print(f"{format_np(multiallelic_DEL_STRs_before_validation_step, DEL_STRs_before_validation_step + multiallelic_DEL_STRs_before_validation_step)} "
          f"are multi-allelic DEL")
    print(f"{format_np(mixed_multiallelic_STRs_before_validation_step, total)} "
          f"mixed INS/DEL multi-allelic STRs in {prefix}")


look_up_liftover_stats("step2:pure_STR:output:", total_STR_variants_before_validation_step)

#%%

# STR Validation Stats - 2nd row of figure: failed liftover

print("\nLiftover failed for:")
liftover_failed_IndelStraddlesMultipleIntevals = int(
    re.search(f"([0-9,]+) IndelStraddlesMultipleIntevals", step1_log_contents).group(1).replace(",", ""))

liftover_failed_for_other_reason = 0
liftover_failed_for_other_reason += int(
    re.search(f"([0-9,]+) CannotLiftOver", step1_log_contents).group(1).replace(",", ""))
liftover_failed_for_other_reason += int(
    re.search(f"([0-9,]+) MismatchedRefAllele", step1_log_contents).group(1).replace(",", ""))
liftover_failed_for_other_reason += int(
    re.search(f"([0-9,]+) NoTarget", step1_log_contents).group(1).replace(",", ""))

print(f"{format_np(liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_np(liftover_failed_for_other_reason, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover for other reasons")

#%%
# STR Validation Stats - 2nd row of figure: liftover succeeded

total_STR_variants_after_1st_liftover_step = int(re.search(
    f"step3:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

print("\nLiftover succeeded for:")
print(f"{format_np(total_STR_variants_after_1st_liftover_step, total_STR_variants_before_validation_step)} "
      f"total STR variants remaining after 1st liftover step")

print("\nLiftover succeeded for (detail):")
look_up_liftover_stats("step3:pure_STR:output:", total_STR_variants_after_1st_liftover_step)


#%%