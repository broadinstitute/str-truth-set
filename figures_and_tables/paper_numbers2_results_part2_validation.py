
import pandas as pd
from figures_and_tables.numbers_utils import format_n, format_np, search

with open("step_A.log", "rt") as f:
    stepA_log_contents = f.read()

#%%

total_STR_variants_before_validation_step = search(
    f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)
total_STR_alleles_before_validation_step = search(
    f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents, type=int)

#%%
print("-"*100)
print("Validation diagram - levels 1, 2:")
print("-"*100)


# STR Validation Stats - 1st row of figure:

prefix = "step2:STR:output:"
total = total_STR_variants_before_validation_step

print(f"{format_n(total_STR_variants_before_validation_step)}    TOTAL STR loci before validation")

INS_STRs_before_validation_step = search(
    f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] INS variants", stepA_log_contents, type=int)

DEL_STRs_before_validation_step = search(
    f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] DEL variants", stepA_log_contents, type=int)

multiallelic_INS_STRs_before_validation_step = search(
    f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] multiallelic INS", stepA_log_contents, type=int)
multiallelic_DEL_STRs_before_validation_step = search(
    f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] multiallelic DEL", stepA_log_contents, type=int)
mixed_multiallelic_STRs_before_validation_step = search(
    f"{prefix}[ ]*([0-9,]+)[ ]* out of[ ]* {total:,d}.*[)] mixed multiallelic", stepA_log_contents, type=int)

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



#%%

# STR Validation Stats - 2nd row of figure: failed liftover

print("\nLiftover failed for:")
liftover_failed_IndelStraddlesMultipleIntevals = search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) IndelStraddlesMultipleIntevals", stepA_log_contents, type=int)

liftover_failed_for_other_reason = 0
liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) CannotLiftOver", stepA_log_contents, type=int)
liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) MismatchedRefAllele", stepA_log_contents, type=int)
liftover_failed_for_other_reason += search(f"step3:STR:failed-liftover:[ ]+([0-9,]+) NoTarget", stepA_log_contents, type=int)

print(f"{format_np(liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_np(liftover_failed_for_other_reason, total_STR_variants_before_validation_step)} "
      f"STRs failed hg38 => T2T liftover for other reasons")

#%%
# STR Validation Stats - 2nd row of figure: liftover succeeded
print("Validation Diagram - level 2:")

total_STR_variants_after_1st_liftover_step = search(
    f"step3:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

print("\nLiftover succeeded for:")
print(f"{format_np(total_STR_variants_after_1st_liftover_step, total_STR_variants_before_validation_step)} "
      f"total STR variants remaining after 1st liftover step")

#print("\nLiftover succeeded for (detail):")
#look_up_liftover_stats("step3:STR:output:", total_STR_variants_after_1st_liftover_step)


#%%
print("-"*100)
print("Validation diagram - level 2 - details:")
print("-"*100)


DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  DEL variants", stepA_log_contents, type=int)
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic DEL variants", stepA_log_contents, type=int)
mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  mixed multiallelic INS/DEL variants", stepA_log_contents, type=int)
assert search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  INS variants", stepA_log_contents, expected_number_of_matches=0) is None
INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0
assert search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic INS variants", stepA_log_contents, expected_number_of_matches=0) is None
multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = 0

print(f"{format_n(DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"DEL variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"MIXED variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")
print(f"{format_n(INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)} "
      f"INS variants failed hg38 to T2T liftover due to IndelStraddlesMultipleIntevals error")


#%%

INS_variants_failed_hg38_to_t2t_liftover = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  INS variants", stepA_log_contents, type=int)
DEL_variants_failed_hg38_to_t2t_liftover = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  DEL variants", stepA_log_contents, type=int)
multiallelic_INS_variants_failed_hg38_to_t2t_liftover = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic INS variants", stepA_log_contents, type=int)
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  multiallelic DEL variants", stepA_log_contents, type=int)
mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filtered:  mixed multiallelic INS/DEL variants", stepA_log_contents, type=int)

total_DEL_variants_failed_hg38_to_t2t_liftover_due_to_other_errors =  DEL_variants_failed_hg38_to_t2t_liftover + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover - DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals
print(f"{format_n(total_DEL_variants_failed_hg38_to_t2t_liftover_due_to_other_errors)} "
      f"DEL variants failed hg38 to T2T liftover due to other errors")
total_mixed_multiallelic_variants_failed_hg38_to_t2t_liftover_due_to_other_errors = mixed_multiallelic_INS_DEL_variants_failed_hg38_to_t2t_liftover - mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals
print(f"{format_n(total_mixed_multiallelic_variants_failed_hg38_to_t2t_liftover_due_to_other_errors)} "
      f"MIXED variants failed hg38 to T2T liftover due to other errors")
total_INS_variants_failed_hg38_to_t2t_liftover_due_to_other_errors = INS_variants_failed_hg38_to_t2t_liftover + multiallelic_INS_variants_failed_hg38_to_t2t_liftover - INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals - multiallelic_INS_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals
print(f"{format_n(total_INS_variants_failed_hg38_to_t2t_liftover_due_to_other_errors)} "
      f"INS variants failed hg38 to T2T liftover due to other errors")

print(F"{format_np(total_DEL_variants_failed_hg38_to_t2t_liftover_due_to_other_errors + total_mixed_multiallelic_variants_failed_hg38_to_t2t_liftover_due_to_other_errors + total_INS_variants_failed_hg38_to_t2t_liftover_due_to_other_errors, total)} "
      F"total variants failed hg38 to T2T liftover due to other errors")
#%%

print("-"*100)
print("Validation diagram - level 3:")
print("-"*100)


# STR Validation Stats - 3rd row of figure: compare alleles with the T2T reference

kept_het_variants = search(
    f"step4:STR[ ]+([0-9,]+) .*kept variants: heterozygous reference genotype", stepA_log_contents, type=int)

kept_variants_insertions_that_match_adjacent_reference_sequence = search(
    f"step4:STR[ ]+([0-9,]+) .*kept variants: insertion matches the adjacent reference sequence", stepA_log_contents, type=int)

total_variants_passed_t2t_comparison = search(
    f"step4:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

assert total_variants_passed_t2t_comparison == kept_het_variants + kept_variants_insertions_that_match_adjacent_reference_sequence

print(f"{format_np(total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
      "Total ")

print(f"{format_np(total_STR_variants_after_1st_liftover_step - total_variants_passed_t2t_comparison, total_STR_variants_before_validation_step)} "
      "failed T2T comparison ")


#%%

print("-"*100)
print("Validation diagram - level 4:")
print("-"*100)

liftover_failed__no_target = search(f"step5:STR:failed-liftover:[ ]+([0-9,]+)[ ]+NoTarget", stepA_log_contents, type=int)

total_variants_passed_liftover_back_to_hg38 = search(f"step5:STR:output:[ ]+([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

print(f"{format_np(liftover_failed__no_target, total_STR_variants_before_validation_step)} "
      "Failed T2T => hg38 liftover: No Target")
print(f"{format_np(total_variants_passed_liftover_back_to_hg38, total_STR_variants_before_validation_step)} "
      "Total passed T2T => hg38 liftover")

total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers = search(
    f"([0-9,]+) .*variants had a different position after hg38 => T2T => hg38", stepA_log_contents,
    expected_number_of_matches=1, type=int)
print(f"{format_n(total_variants_with_different_position_after_hg38_to_t2t_to_hg38_liftovers)} "
      f"variants had a different position after hg38 => T2T => hg38 liftovers")

#%%
print("-"*100)
print("Validation diagram - final level:")
print("-"*100)

total_variants_passed_all_validation_steps = search(
    f"step7:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)

print(f"{format_np(total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants passed all validation steps")

print(f"{format_np(total_STR_variants_before_validation_step - total_variants_passed_all_validation_steps, total_STR_variants_before_validation_step)} "
      "Total variants failed some validation step")

#%%
print("-"*100)
print("Validation diagram - final level - details:")
print("-"*100)


df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")

INS_variants_in_final_truth_set = search(
    f"step7:STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] INS variants", stepA_log_contents, type=int)
DEL_variants_in_final_truth_set = search(
    f"step7:STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] DEL variants", stepA_log_contents, type=int)
multiallelic_INS_variants_in_final_truth_set = search(
    f"step7:STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic INS variants", stepA_log_contents, type=int)
multiallelic_DEL_variants_in_final_truth_set = search(
    f"step7:STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] multiallelic DEL variants", stepA_log_contents, type=int)
mixed_multiallelic_INS_DEL_variants_in_final_truth_set = search(
    f"step7:STR:output:[ ]*([0-9,]+) out of [ ]*{len(df_variants):,d}.*[%][)] mixed multiallelic INS/DEL variants", stepA_log_contents, type=int)

print(f"{format_n(DEL_variants_in_final_truth_set + multiallelic_DEL_variants_in_final_truth_set)} "
      f"DEL variants in final truth set")
print(f"{format_n(mixed_multiallelic_INS_DEL_variants_in_final_truth_set)} "
      f"MIXED variants in final truth set")
print(f"{format_n(INS_variants_in_final_truth_set + multiallelic_INS_variants_in_final_truth_set)} "
      f"INS variants in final truth set")

#%%
print("-"*100)
print("Validation diagram - overall percent that could be validated and that were validated:")
print("-"*100)
total_STR_variants_that_could_be_validated = total_STR_variants_before_validation_step - liftover_failed_IndelStraddlesMultipleIntevals
print(f"{format_n(total_STR_variants_that_could_be_validated)} STRs that could be validated against T2T")
print(f"{format_np(len(df_variants) - liftover_failed_IndelStraddlesMultipleIntevals, total_STR_variants_that_could_be_validated)} "
      f"STRs passed this validation")


#%%


#df_variants_before_validation = pd.read_table("step2.STRs.variants.tsv.gz")
#df_insertions = df_variants_before_validation[df_variants_before_validation.INS_or_DEL.isin({"INS", "INS:INS"})]
print(f"{INS_STRs_before_validation_step + multiallelic_INS_STRs_before_validation_step:,d} monoallellic and multiallelic expansions before validation")

print(f"{format_np(INS_variants_in_final_truth_set + multiallelic_INS_variants_in_final_truth_set, INS_STRs_before_validation_step + multiallelic_INS_STRs_before_validation_step)} "
      f"monoallellic and multiallelic expansions passed validation")

print(f"liftover didn't work for {format_np(DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals+multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals, DEL_STRs_before_validation_step)} contractions ")
print(f"liftover didn't work for {format_np(mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals,mixed_multiallelic_STRs_before_validation_step)} mixed multiallelics")

liftover_did_work_for = mixed_multiallelic_STRs_before_validation_step + DEL_STRs_before_validation_step+multiallelic_DEL_STRs_before_validation_step - liftover_failed_IndelStraddlesMultipleIntevals
n_passed = (DEL_variants_in_final_truth_set + mixed_multiallelic_INS_DEL_variants_in_final_truth_set
            - mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals
            - DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals)
print(f"Liftover did work for {format_n(liftover_did_work_for)} contraction and mixed multiallelics.")
print(f"Of these {format_np(n_passed, liftover_did_work_for)} passed all validation steps.")
#%%