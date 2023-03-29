
import collections
import gzip
import pandas as pd
import numpy as np
from figures_and_tables.numbers_utils import format_n, format_np, search



#%%
print("---")
print("Numbers from the Abstract:")

with open("step_A.log", "rt") as f:
    stepA_log_contents = f.read()


total_variants_in_syndip = search(f"[ ]+([0-9,]+)[ ]+total[ ]variants",
                                  stepA_log_contents, expected_number_of_matches=2, use_match_i=0, type=int)
total_indels_in_syndip = search(f"[ ]+([0-9,]+)[ ]+INDELs",
                                stepA_log_contents, expected_number_of_matches=4, use_match_i=0, type=int)
high_confidence_variants_in_syndip = search(f"[ ]+([0-9,]+)[ ]+total[ ]variants",
                                            stepA_log_contents, expected_number_of_matches=2, use_match_i=1, type=int)
high_confidence_indels_in_syndip = search(f"[ ]+([0-9,]+)[ ]+INDELs",
                                          stepA_log_contents, expected_number_of_matches=4, use_match_i=2, type=int)

print(f"{format_n(total_variants_in_syndip)} total variants in SynDip")
print(f"{format_np(total_indels_in_syndip, total_variants_in_syndip)} total INDELs in SynDip")
print("")
print(f"{format_n(high_confidence_variants_in_syndip)} high-confidence variants in SynDip")
print(f"{format_np(high_confidence_indels_in_syndip, high_confidence_variants_in_syndip)} high-confidence INDELs in SynDip")


#%%
print("---")
all_high_confidence_regions = []
with gzip.open("./ref/full.38.bed.gz", "rt") as f:
    for line in f:
        fields = line.strip().split("\t")
        all_high_confidence_regions.append(int(fields[2]) - int(fields[1]))


print(f"{format_n(sum(all_high_confidence_regions))}bp total size of high confidence regions")
print(f"{format_n(int(np.mean(all_high_confidence_regions)))}bp mean size of high confidence regions")
print(f"{format_n(len(all_high_confidence_regions))} high confidence regions")

#%%
print("---")
df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_pure_variants = df_variants[df_variants.IsPureRepeat]
df_interrupted_variants = df_variants[~df_variants.IsPureRepeat]

total_variants = len(df_variants)
print(f"{format_n(total_variants)} total variants in the TR truth set")
print(f"{format_np(len(df_pure_variants), total_variants)} pure repeats")
print(f"{format_np(len(df_interrupted_variants), total_variants)} repeats with interruptions")


#%%
print("---")
total_TR_variants_before_validation_step = search(f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL variants",
                                                   stepA_log_contents, type=int)
total_TR_alleles_before_validation_step = search(f"step2:STR:output:[ ]*([0-9,]+)[ ]* TOTAL alleles",
                                                  stepA_log_contents, type=int)

print(f"{total_TR_variants_before_validation_step:10,d} total TR variants before validation")
print(f"{100*total_TR_variants_before_validation_step/high_confidence_indels_in_syndip:9.0f}% "
      f"total TR variants before validation as % of high-confidence indels")

#%%

df_variants_before_validation = pd.read_table("step2.STRs.variants.tsv.gz")
df_variants_before_validation_3_to_24bp_motifs = df_variants_before_validation[
    (df_variants_before_validation.MotifSize >= 3) & (df_variants_before_validation.MotifSize <= 24)
]

total = len(df_variants_before_validation)
x = sum(df_variants_before_validation.MotifSize == 2)
print(f"{x:,d} TR loci ({100*x/total:,.0f}%) have motif sizes == 2bp")
y = len(df_variants_before_validation_3_to_24bp_motifs)
print(f"{y:,d} TR loci ({100*y/total:,.0f}%) have motif sizes >= 3bp and <= 24bp")
z = sum(df_variants_before_validation.MotifSize > 24)
print(f"{z:,d} TR loci ({100*z/total:,.0f}%) have motif sizes > 24bp")

total2 = len(df_variants_before_validation_3_to_24bp_motifs)
a = sum(df_variants_before_validation_3_to_24bp_motifs.IsPureRepeat)
print(f"{a:,d} TR loci ({100*a/total2:,.0f}%) of loci with 3-24bp motifs are pure repeat")
b = sum(~df_variants_before_validation_3_to_24bp_motifs.IsPureRepeat)
print(f"{b:,d} TR loci ({100*b/total2:,.0f}%) of loci with 3-24bp motifs have interruptions")


#%%
print("---")
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_pure_alleles = df_alleles[df_alleles.IsPureRepeat]
df_interrupted_alleles = df_alleles[~df_alleles.IsPureRepeat]
print(f"{len(df_alleles):,d} TR alleles")
print(f"{len(df_pure_alleles):,d} pure alleles")
print(f"{len(df_interrupted_alleles):,d} interrupted alleles")


#%%