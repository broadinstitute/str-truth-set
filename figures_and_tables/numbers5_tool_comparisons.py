import pandas as pd
from figures_and_tables.numbers_utils import format_np

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat == "Yes"]
df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_alleles = df_alleles[df_alleles.IsPureRepeat == "Yes"]

df_alleles_tool_comparison = pd.read_table("./tool_comparison/combined.results.alleles.tsv.gz")
df_alleles_tool_comparison = df_alleles_tool_comparison[df_alleles_tool_comparison.IsPureRepeat]

print("\n".join(df_alleles_tool_comparison.columns))

#%%

print("% STR truth set loci with 2-6bp motifs")
df_variants_filtered = df_variants[(df_variants.MotifSize <= 6) & (df_variants.IsFoundInReference == "Yes")]

print(format_np(len(df_variants_filtered), len(df_variants)), "STR truth set loci with 2-6bp motifs")

#%%

print("Exact accuracy: ")

print("    -----")
df_alleles_tool_comparison_filtered = df_alleles_tool_comparison[
    (df_alleles_tool_comparison.IsFoundInReference) &
    (df_alleles_tool_comparison.PositiveOrNegative == "positive") &
    (df_alleles_tool_comparison.MotifSize <= 6) &
    (df_alleles_tool_comparison.MotifSize >= 2)
]

for tool in ["ExpansionHunter", "GangSTR", "HipSTR"]:
    for coverage in "40x", "20x", "10x", "05x":
        df_alleles_tool_comparison_filtered_by_coverage = df_alleles_tool_comparison_filtered[
            (df_alleles_tool_comparison.coverage == coverage) &
            (df_alleles_tool_comparison["DiffFromRefRepeats: Allele: Truth"] != 0)
        ]
        column = f"DiffRepeats: Allele: {tool} - Truth"
        exactly_right_genotypes_count = sum(df_alleles_tool_comparison_filtered_by_coverage[column] == 0)
        df_na = df_alleles_tool_comparison_filtered_by_coverage[df_alleles_tool_comparison_filtered_by_coverage[column].isna()]
        df_not_na = df_alleles_tool_comparison_filtered_by_coverage[df_alleles_tool_comparison_filtered_by_coverage[column].notna()]
        print(format_np(exactly_right_genotypes_count, len(df_alleles_tool_comparison_filtered_by_coverage)), f"of STR "
              f"calls by {tool} are exactly right at {coverage} coverage, but",
              format_np(len(df_na), len(df_alleles_tool_comparison_filtered_by_coverage)), " alleles are NA",
              format_np(len(set(df_not_na.LocusId)), len(set(df_alleles_tool_comparison_filtered_by_coverage.LocusId))), " loci are not NA",
        )
    print("---")

#%%
df_alleles_tool_comparison_filtered_by_coverage = df_alleles_tool_comparison_filtered[
    (df_alleles_tool_comparison.coverage == "40x") &
    (df_alleles_tool_comparison["DiffFromRefRepeats: Allele: Truth"] != 0)
]
print("% STR truth set loci not called by HipSTR")

print(len(set(df_alleles_tool_comparison_filtered_by_coverage.LocusId)))
len(set(df_alleles_tool_comparison_filtered_by_coverage[df_alleles_tool_comparison_filtered_by_coverage[f"DiffRepeats: Allele: HipSTR - Truth"].isna()].LocusId))
#collections.Counter(len(set([1])))

#%%
sum(df_alleles_tool_comparison_filtered_by_coverage[f"DiffRepeats: Allele: HipSTR - Truth"].isna() & (df_alleles_tool_comparison_filtered_by_coverage[f"DiffRepeats: Allele: HipSTR - Truth"] != 0))


#%%