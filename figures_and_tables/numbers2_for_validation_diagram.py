"""This script computes numbers that appear in the truth set validation flow chart"""

import pandas as pd
from figures_and_tables.numbers_utils import format_n, format_np, search

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat == "Yes"]

#%%

with open("step_A.log", "rt") as f:
    stepA_log_contents = f.read()

#%%

total_STR_variants_before_validation_step = int(search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))
total_STR_alleles_before_validation_step = int(search(
    f"step2:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents).replace(",", ""))

#%%
print("-"*100)
print("Validation diagram - levels 1, 2:")
print("-"*100)


# STR Validation Stats - 1st row of figure:

print(f"{format_n(total_STR_variants_before_validation_step)}    TOTAL STR loci before validation")

def look_up_liftover_stats(prefix, total):
    INS_STRs_before_validation_step = int(search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] INS variants", stepA_log_contents).replace(",", ""))

    DEL_STRs_before_validation_step = int(search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] DEL variants", stepA_log_contents).replace(",", ""))

    multiallelic_INS_STRs_before_validation_step = int(search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] multiallelic INS", stepA_log_contents).replace(",", ""))
    multiallelic_DEL_STRs_before_validation_step = int(search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] multiallelic DEL", stepA_log_contents).replace(",", ""))
    mixed_multiallelic_STRs_before_validation_step = int(search(
        f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] mixed multiallelic", stepA_log_contents).replace(",", ""))

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
    search(f"([0-9,]+) IndelStraddlesMultipleIntevals", stepA_log_contents, expected_number_of_matches=2).replace(",", ""))

liftover_failed_for_other_reason = 0
liftover_failed_for_other_reason += int(
    search(f"([0-9,]+) CannotLiftOver", stepA_log_contents, expected_number_of_matches=2).replace(",", ""))
liftover_failed_for_other_reason += int(
    search(f"([0-9,]+) MismatchedRefAllele", stepA_log_contents, expected_number_of_matches=2).replace(",", ""))
liftover_failed_for_other_reason += int(
    search(f"([0-9,]+) NoTarget", stepA_log_contents, expected_number_of_matches=2).replace(",", ""))

print(f"{format_np(liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_np(liftover_failed_for_other_reason, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover for other reasons")

#%%
# STR Validation Stats - 2nd row of figure: liftover succeeded
print("Validation Diagram - level 2:")

total_STR_variants_after_1st_liftover_step = int(search(
    f"step3:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))

print("\nLiftover succeeded for:")
print(f"{format_np(total_STR_variants_after_1st_liftover_step, total_STR_variants_before_validation_step)} "
      f"total STR variants remaining after 1st liftover step")

#print("\nLiftover succeeded for (detail):")
#look_up_liftover_stats("step3:pure_STR:output:", total_STR_variants_after_1st_liftover_step)


#%%
print("-"*100)
print("Validation diagram - level 2 - details:")
print("-"*100)


DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  DEL variants", stepA_log_contents).replace(",", ""))
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic DEL variants", stepA_log_contents).replace(",", ""))
mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  mixed multiallelic INS/DEL variants", stepA_log_contents).replace(",", ""))
assert search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  INS variants", stepA_log_contents, expected_number_of_matches=0) is None
INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0
assert search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic INS variants", stepA_log_contents, expected_number_of_matches=0) is None
multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0

print(f"{format_n(DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"DEL variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"MIXED variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"INS variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")


#%%

INS_variants_failed_hg38_to_t2t_liftover = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  INS variants", stepA_log_contents).replace(",", ""))
DEL_variants_failed_hg38_to_t2t_liftover = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  DEL variants", stepA_log_contents).replace(",", ""))
multiallelic_INS_variants_failed_hg38_to_t2t_liftover = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic INS variants", stepA_log_contents).replace(",", ""))
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic DEL variants", stepA_log_contents).replace(",", ""))
mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover = int(search(
    f"step3:pure_STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  mixed multiallelic INS/DEL variants", stepA_log_contents).replace(",", ""))

print(f"{format_n(DEL_variants_failed_hg38_to_t2t_liftover + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover - DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"DEL variants failed hg38 to T2T liftover due to other errors")
print(f"{format_n(mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover - mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"MIXED variants failed hg38 to T2T liftover due to other errors")
print(f"{format_n(INS_variants_failed_hg38_to_t2t_liftover + multiallelic_INS_variants_failed_hg38_to_t2t_liftover - INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"INS variants failed hg38 to T2T liftover due to other errors")

#%%

print("-"*100)
print("Validation diagram - level 3:")
print("-"*100)


# STR Validation Stats - 3rd row of figure: compare alleles with the T2T reference

kept_het_variants = int(search(
    f"([0-9,]+) .*kept variants: heterozygous reference genotype", stepA_log_contents).replace(",", ""))

kept_variants_insertions_that_match_adjacent_reference_sequence = int(search(
    f"([0-9,]+) .*kept variants: insertion matches the adjacent reference sequence", stepA_log_contents).replace(",", ""))

total_variants_passed_t2t_comparison = int(search(
    f"step4:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))

assert total_variants_passed_t2t_comparison == kept_het_variants + kept_variants_insertions_that_match_adjacent_reference_sequence

print(f"{format_np(total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
    "Total ")

print(f"{format_np(total_STR_variants_after_1st_liftover_step - total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
      "failed T2T comparison ")


#%%

print("-"*100)
print("Validation diagram - level 4:")
print("-"*100)

liftover_failed__no_target = int(
    search(f"step5:pure_STR:failed-liftover:[ ]*([0-9,]+).*NoTarget", stepA_log_contents).replace(",", ""))

total_variants_passed_liftover_back_to_hg38 = int(search(
    f"step5:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))

print(f"{format_np(liftover_failed__no_target, total_STR_variants_before_validation_step)} "
      "Failed T2T => hg38 liftover: No Target")
print(f"{format_np(total_variants_passed_liftover_back_to_hg38, total_STR_variants_before_validation_step)} "
      "Total passed T2T => hg38 liftover")

total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers = int(search(
    f"([0-9,]+) .*variants had a different position after hg38 => T2T => hg38", stepA_log_contents).replace(",", ""))
print(f"{format_n(total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers)} "
      f"variants had a different position after hg38 => T2T => hg38 liftovers")

#%%
print("-"*100)
print("Validation diagram - final level:")
print("-"*100)

total_variants_passed_all_validation_steps = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents).replace(",", ""))

print(f"{format_np(total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants passed all validation steps")

print(f"{format_np(total_STR_variants_before_validation_step - total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants failed some validation step")

#%%
print("-"*100)
print("Validation diagram - final level - details:")
print("-"*100)

INS_variants_in_final_truth_set = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] INS variants", stepA_log_contents).replace(",", ""))
DEL_variants_in_final_truth_set = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] DEL variants", stepA_log_contents).replace(",", ""))
multiallelic_INS_variants_in_final_truth_set = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic INS variants", stepA_log_contents).replace(",", ""))
multiallelic_DEL_variants_in_final_truth_set = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic DEL variants", stepA_log_contents).replace(",", ""))
mixed_multiallelic_INS_DEL_variants_in_final_truth_set = int(search(
    f"step7:pure_STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] mixed multiallelic INS/DEL variants", stepA_log_contents).replace(",", ""))

print(f"{format_n(DEL_variants_in_final_truth_set + multiallelic_DEL_variants_in_final_truth_set)} "
      f"DEL variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(mixed_multiallelic_INS_DEL_variants_in_final_truth_set)} "
      f"MIXED variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(INS_variants_in_final_truth_set + multiallelic_INS_variants_in_final_truth_set)} "
      f"INS variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")

#%%
print("-"*100)
print("Validation diagram - overall percent that could be validated and that were validated:")
print("-"*100)
total_STR_variants_that_could_be_validated = total_STR_variants_before_validation_step - liftover_failed_IndelStraddlesMultipleIntevals
print(f"{format_n(total_STR_variants_that_could_be_validated)} STRs that could be validated against T2T")
print(f"{format_np(len(df_variants) - liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_that_could_be_validated)} "
      f"STRs passed this validation")


#%%