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
print("-"*100)
print("Validation vs. T2T:")
print("-"*100)


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

#print("\nLiftover succeeded for (detail):")
#look_up_liftover_stats("step3:pure_STR:output:", total_STR_variants_after_1st_liftover_step)


#%%

# STR Validation Stats - 3rd row of figure: compare alleles with the T2T reference

kept_het_variants = int(re.search(
    f"([0-9,]+) .*kept variants: heterozygous reference genotype", step1_log_contents).group(1).replace(",", ""))

kept_variants_insertions_that_match_adjacent_reference_sequence = int(re.search(
    f"([0-9,]+) .*kept variants: insertion matches the adjacent reference sequence", step1_log_contents).group(1).replace(",", ""))

total_variants_passed_t2t_comparison = int(re.search(
    f"step4:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

assert total_variants_passed_t2t_comparison == kept_het_variants + kept_variants_insertions_that_match_adjacent_reference_sequence

print(f"{format_np(total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
    "Total ")

print(f"{format_np(total_STR_variants_after_1st_liftover_step - total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
      "failed t2t comparison ")


#%%

total_variants_passed_liftover_back_to_hg38 = int(re.search(
    f"step5:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

print(f"{format_np(total_variants_passed_liftover_back_to_hg38, total_STR_variants_before_validation_step)} "
      "Total passed t2t => hg38 liftover")

print("-"*100)

total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers = int(re.search(
    f"([0-9,]+) .*variants had a different position after hg38 => T2T => hg38", step1_log_contents).group(1).replace(",", ""))
print(total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers, "variants had a different position after hg38 => T2T => hg38")

#%%
print("-"*100)

total_variants_passed_all_validation_steps = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", step1_log_contents).group(1).replace(",", ""))

print(f"{format_np(total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants passed all validation steps")

print(f"{format_np(total_STR_variants_before_validation_step - total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants failed some validation step")

#%%
print("-"*100)

DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  DEL variants", step1_log_contents).group(1).replace(",", ""))
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic DEL variants", step1_log_contents).group(1).replace(",", ""))
mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  mixed multiallelic INS/DEL variants", step1_log_contents).group(1).replace(",", ""))
assert re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  INS variants", step1_log_contents) is None
INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0
assert re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic INS variants", step1_log_contents) is None
multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0

print(f"{format_n(DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"DEL variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"MIXED variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"INS variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")


#%%

print("-"*100)

INS_variants_failed_hg38_to_t2t_liftover = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  INS variants", step1_log_contents).group(1).replace(",", ""))
DEL_variants_failed_hg38_to_t2t_liftover = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  DEL variants", step1_log_contents).group(1).replace(",", ""))
multiallelic_INS_variants_failed_hg38_to_t2t_liftover = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic INS variants", step1_log_contents).group(1).replace(",", ""))
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic DEL variants", step1_log_contents).group(1).replace(",", ""))
mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover = int(re.search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  mixed multiallelic INS/DEL variants", step1_log_contents).group(1).replace(",", ""))

print(f"{format_n(DEL_variants_failed_hg38_to_t2t_liftover + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover - DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"DEL variants failed hg38 to T2T liftover due to other errors")
print(f"{format_n(mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover - mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"MIXED variants failed hg38 to T2T liftover due to other errors")
print(f"{format_n(INS_variants_failed_hg38_to_t2t_liftover + multiallelic_INS_variants_failed_hg38_to_t2t_liftover - INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"INS variants failed hg38 to T2T liftover due to other errors")

#%%

INS_variants_in_final_truthset = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] INS variants", step1_log_contents).group(1).replace(",", ""))
DEL_variants_in_final_truthset = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] DEL variants", step1_log_contents).group(1).replace(",", ""))
multiallelic_INS_variants_in_final_truthset = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic INS variants", step1_log_contents).group(1).replace(",", ""))
multiallelic_DEL_variants_in_final_truthset = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic DEL variants", step1_log_contents).group(1).replace(",", ""))
mixed_multiallelic_INS_DEL_variants_in_final_truthset = int(re.search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] mixed multiallelic INS/DEL variants", step1_log_contents).group(1).replace(",", ""))

print(f"{format_n(DEL_variants_in_final_truthset + multiallelic_DEL_variants_in_final_truthset)} "
      f"DEL variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(mixed_multiallelic_INS_DEL_variants_in_final_truthset)} "
      f"MIXED variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(INS_variants_in_final_truthset + multiallelic_INS_variants_in_final_truthset)} "
      f"INS variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")

