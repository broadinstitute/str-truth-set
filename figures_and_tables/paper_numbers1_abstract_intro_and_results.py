
import gzip
import pandas as pd
import numpy as np
from figures_and_tables.numbers_utils import format_n, format_np, search

print("---")
print("Numbers from the Abstract: ")

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_pure_variants = df_variants[df_variants.IsPureRepeat]
df_interrupted_variants = df_variants[~df_variants.IsPureRepeat]

total_variants = len(df_variants)
print(f"{format_n(total_variants)} total variants in the TR truth set")
print(f"{format_np(len(df_pure_variants), total_variants)} pure variants")
print(f"{format_np(len(df_interrupted_variants), total_variants)} variants with interruptions")


#%%
print("---")
print("Numbers from Results Section 1: Deriving a TR truth set from the Synthetic Diploid Benchmark")

all_high_confidence_regions = []
with gzip.open("./ref/full.38.bed.gz", "rt") as f:
    for line in f:
        fields = line.strip().split("\t")
        all_high_confidence_regions.append(int(fields[2]) - int(fields[1]))

chrY_size = 57_227_415   # based on https://www.ncbi.nlm.nih.gov/grc/human/data
hg38_genome_total_size = 3_088_269_832 - chrY_size

print(f"{format_np(sum(all_high_confidence_regions), hg38_genome_total_size)} bp total size of high confidence regions")
#print(f"{format_n(int(np.mean(all_high_confidence_regions)))}bp mean size of high confidence regions")
#print(f"{format_n(len(all_high_confidence_regions))} high confidence regions")
print("")

with open("step_A.log", "rt") as f:
    stepA_log_contents = f.read()

total_variants_in_syndip = search(f"step1:input:[ ]+([0-9,]+)[ ]+TOTAL[ ]variants", stepA_log_contents, type=int)
total_alleles_in_syndip = search(f"step1:input:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)

high_confidence_variants_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]variants", stepA_log_contents, type=int)
high_confidence_alleles_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)

high_confidence_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] SNV variants", stepA_log_contents).replace(",", ""))
high_confidence_multiallelic_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] multiallelic SNV variants", stepA_log_contents).replace(",", ""))
high_confidence_indel_variants_in_syndip = high_confidence_variants_in_syndip - high_confidence_SNV_variants - high_confidence_multiallelic_SNV_variants

high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))
high_confidence_indel_alleles_in_syndip = high_confidence_INS_alleles + high_confidence_DEL_alleles

print(f"{format_n(total_variants_in_syndip)} total variants in SynDip")
#print(f"{format_n(total_alleles_in_syndip)} total alleles in SynDip")
print(f"{format_n(high_confidence_variants_in_syndip)} high-confidence variants in SynDip")
#print(f"{format_n(high_confidence_alleles_in_syndip)} high-confidence alleles in SynDip")

print(f"{format_np(high_confidence_indel_variants_in_syndip, high_confidence_variants_in_syndip)} high-confidence INDEL variants in SynDip")
#print(f"{format_np(high_confidence_indel_alleles_in_syndip, high_confidence_alleles_in_syndip)} high-confidence INDEL alleles in SynDip")



#%%

print("---")
print("Numbers from Results Section 2: Validating the TR truth set using the telomere-to-telomere reference genome")

total_TR_variants_before_validation_step = search(f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants", stepA_log_contents, type=int)
total_TR_alleles_before_validation_step = search(f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles", stepA_log_contents, type=int)

#print(f"{total_TR_variants_before_validation_step:10,d} total TR variants before validation")
#print(f"{100*total_TR_variants_before_validation_step/high_confidence_indel_variants_in_syndip:9.0f}% "
#      f"total TR variants before validation as % of high-confidence indels")

df_variants_before_validation = pd.read_table("step2.STRs.variants.tsv.gz")
#df_alleles_before_validation = pd.read_table("step2.STRs.alleles.tsv.gz")

if total_TR_variants_before_validation_step != len(df_variants_before_validation):
    raise ValueError(f"total_TR_variants_before_validation_step != len(df_variants_before_validation)")

expansion_variants_before_validation = sum(df_variants_before_validation.INS_or_DEL.isin(('INS', 'INS:INS')))
contraction_variants_before_validation = sum(df_variants_before_validation.INS_or_DEL.isin(('DEL', 'DEL:DEL')))
mixed_variants_before_validation = sum(df_variants_before_validation.INS_or_DEL == 'DEL:INS')
print(f"{format_np(expansion_variants_before_validation, len(df_variants_before_validation))} TR expansion variants before validation")
print(f"{format_np(contraction_variants_before_validation, len(df_variants_before_validation))} TR contraction variants before validation")
print(f"{format_np(mixed_variants_before_validation, len(df_variants_before_validation))} TR mixed variants before validation")

#print(f"{format_np(len(df_alleles_before_validation), high_confidence_alleles_in_syndip)} TR alleles before validation")

expansion_variants_after_validation = sum(df_variants.INS_or_DEL.isin(('INS', 'INS:INS')))
contraction_variants_after_validation = sum(df_variants.INS_or_DEL.isin(('DEL', 'DEL:DEL')))
mixed_variants_after_validation = sum(df_variants.INS_or_DEL == 'DEL:INS')
print("-")
print(f"{format_np(expansion_variants_after_validation, expansion_variants_before_validation)} TR expansion variants after validation")
print(f"{format_np(contraction_variants_after_validation, contraction_variants_before_validation)} TR contraction variants after validation")
print(f"{format_np(mixed_variants_after_validation, mixed_variants_before_validation)} TR mixed variants after validation")
print("-")
print(f"{format_np(len(df_variants_before_validation) - len(df_variants), len(df_variants_before_validation))} "
      f"total variants didn't pass validation")




#%%
DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  DEL variants", stepA_log_contents, type=int)
multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  multiallelic DEL variants", stepA_log_contents, type=int)
mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals = search(
    f"step3:STR:failed-liftover:[ ]*([0-9,]+) out of .* filter:IndelStraddlesMultipleIntevals  mixed multiallelic INS/DEL variants", stepA_log_contents, type=int)


print(f"{format_np(DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals + multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals, contraction_variants_before_validation)} "
      f"contraction variants where liftover failed due to IndelStraddlesMultipleIntevals")

print(f"{format_np(mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals, mixed_variants_before_validation)} "
      f"mixed multi-allelic variants where liftover failed due to IndelStraddlesMultipleIntevals")

#%%
print("--")

contraction_and_mixed_variants_with_IndelStraddlesMultipleIntervals = (
        DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals +
        multiallelic_DEL_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals +
        mixed_variants_failed_hg38_to_t2t_liftover_due_to_IndelStraddlesMultipleIntevals
)
contraction_and_mixed_variants_without_IndelStraddlesMultipleIntevals = (
        contraction_variants_before_validation + mixed_variants_before_validation - contraction_and_mixed_variants_with_IndelStraddlesMultipleIntervals
)

contraction_and_mixed_variants_that_passed_validation = (
    contraction_variants_after_validation + mixed_variants_after_validation - contraction_and_mixed_variants_with_IndelStraddlesMultipleIntervals
)

print(f"{format_np(contraction_and_mixed_variants_that_passed_validation, contraction_and_mixed_variants_without_IndelStraddlesMultipleIntevals)} "
      f"contraction and mixed multi-allelic variants that passed validation")

print("--")

print(f"{format_np(expansion_variants_after_validation, len(df_variants))} TR expansion variants after validation")
print(f"{format_np(contraction_variants_after_validation, len(df_variants))} TR contraction variants after validation")
print(f"{format_np(mixed_variants_after_validation, len(df_variants))} TR mixed variants after validation")



#%%


print("Numbers from Results Section 3: Truth set summary")
total_variants = len(df_variants)
print(f"{format_n(high_confidence_indel_variants_in_syndip)} indel variants in SynDip within high-confidence regions")
print(f"{format_np(total_variants, high_confidence_indel_variants_in_syndip)} total variants in the TR truth set")
print(f"{format_np(len(df_pure_variants), total_variants)} pure variants")
print(f"{format_np(len(df_interrupted_variants), total_variants)} variants with interruptions")

print(f"{format_np(sum(~df_variants.IsMultiallelic), total_variants)} mono-allelic variants")
print(f"{format_np(sum(df_variants.IsMultiallelic), total_variants)} multi-allelic variants")

#%%
print("--")
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_alleles["ReferenceAlleleSize (bp)"] = df_alleles.NumRepeatsInReference * df_alleles.MotifSize
df_alleles["AlleleSizeMinusReference (bp)"] = df_alleles["RepeatSize (bp)"] - df_alleles["ReferenceAlleleSize (bp)"]
df_alleles["AlleleSizeMinusReference.abs() (bp)"] = df_alleles["AlleleSizeMinusReference (bp)"].abs()

print(f"{format_np(sum(df_alleles['AlleleSizeMinusReference (bp)'] < -50), len(df_alleles))} contraction alleles > 50bp away from reference")
print(f"{format_np(sum(df_alleles['AlleleSizeMinusReference (bp)'] > 50), len(df_alleles))} expansion alleles > 50bp away from reference")
print(f"{format_np(sum(df_alleles['AlleleSizeMinusReference.abs() (bp)'] > 50), len(df_alleles))} alleles > 50bp away from reference")
print(f"{format_np(sum(df_alleles['AlleleSizeMinusReference.abs() (bp)'] <= 30), len(df_alleles))} alleles <= 30bp")

#%%

print("--")

#['Chrom', 'Start1Based', 'End1Based', 'Locus', 'LocusId', 'INS_or_DEL', 'HET_or_HOM', 'Motif', 'CanonicalMotif', 'MotifSize', 'NumRepeatsInReference', 'VcfPos', 'VcfRef', 'VcfAlt', 'VcfGenotype', 'SummaryString', 'IsFoundInReference', 'IsPureRepeat', 'IsMultiallelic', 'NumRepeats', 'RepeatSize (bp)', 'NumPureRepeats', 'PureRepeatSize (bp)', 'FractionPureRepeats', 'RepeatUnitInterruptionIndex', 'OverlapsIlluminaSTRCatalog: Locus', 'OverlapsIlluminaSTRCatalog: Motif', 'OverlapsGangSTRCatalog17: Locus', 'OverlapsGangSTRCatalog17: Motif', 'OverlapsGangSTRCatalog13: Locus', 'OverlapsGangSTRCatalog13: Motif', 'OverlapsHipSTRCatalog: Locus', 'OverlapsHipSTRCatalog: Motif', 'OverlapsKnownDiseaseAssociatedSTRs: Locus', 'OverlapsKnownDiseaseAssociatedSTRs: Motif', 'OverlapsTRFPureRepeats15bp: Locus', 'OverlapsTRFPureRepeats15bp: Motif', 'OverlapsTRFPureRepeats12bp: Locus', 'OverlapsTRFPureRepeats12bp: Motif', 'OverlapsTRFPureRepeats9bp: Locus', 'OverlapsTRFPureRepeats9bp: Motif', 'OverlapsTRFPureRepeats6bp: Locus', 'OverlapsTRFPureRepeats6bp: Motif', 'OverlapsSegDupIntervals', 'GeneRegionFromGencode_V42', 'GeneNameFromGencode_V42', 'GeneIdFromGencode_V42', 'TranscriptIdFromGencode_V42', 'GeneRegionFromMane_V1', 'GeneNameFromMane_V1', 'GeneIdFromMane_V1', 'TranscriptIdFromMane_V1']

# compute median allele size in hg38 and in the truth set
median_repeats_in_ref_alleles = int(df_alleles.NumRepeatsInReference.median())
median_repeats_in_non_ref_alleles = int(df_alleles.NumRepeats.median())
print(f"{format_n(median_repeats_in_ref_alleles)} median_repeats_in_ref_alleles")
print(f"{format_n(median_repeats_in_non_ref_alleles)} median_repeats_in_non_ref_alleles")

median_reference_allele_size = int(df_alleles["ReferenceAlleleSize (bp)"].median())
median_non_reference_allele_size = int(df_alleles["RepeatSize (bp)"].median())

print(f"{format_n(median_reference_allele_size)} median_reference_allele_size")
print(f"{format_n(median_non_reference_allele_size)} median_non_reference_allele_size")

#%%
