import os
import pandas as pd

tool_comparison_base_dir = "tool_comparison"


def locus_id_to_region(s):
    chrom, start_0based, end, motif = s.split("-")
    return f"{chrom}:{int(start_0based) + 1}-{end}"


variants_tables = []
alleles_tables = []
ehdn_tables = []
for coverage_label, results_directory in [
    ("40x", "results"),
    ("30x", "results_for_downsampled_30x_bam"),
    ("20x", "results_for_downsampled_20x_bam"),
    ("10x", "results_for_downsampled_10x_bam"),
    ("05x", "results_for_downsampled_5x_bam"),
    ("exome", "results_for_exome"),
]:
    # ExpansionHunter, GangSTR, HipSTR
    print(f"Processing {results_directory}")
    for table_list, filename_suffix in [
        (variants_tables, ".variants"),
        (alleles_tables, ".alleles"),
    ]:
        positive_df_path = os.path.join(
            tool_comparison_base_dir, results_directory, f"STR_truth_set.v1.for_comparison{filename_suffix}.tsv")
        print(f"Loading {positive_df_path}")
        positive_df = pd.read_table(positive_df_path, dtype=str)

        negative_df_path = os.path.join(
            tool_comparison_base_dir, results_directory, f"negative_loci.for_comparison{filename_suffix}.tsv")
        print(f"Loading {negative_df_path}")
        negative_df = pd.read_table(negative_df_path, dtype=str)

        positive_df.loc[:, "Locus"] = positive_df.LocusId.apply(locus_id_to_region)
        negative_df.loc[:, "Locus"] = negative_df.LocusId.apply(locus_id_to_region)
        positive_df.loc[:, "Coverage"] = coverage_label
        negative_df.loc[:, "Coverage"] = coverage_label
        positive_df.loc[:, "PositiveOrNegative"] = "positive"
        negative_df.loc[:, "PositiveOrNegative"] = "negative"

        table_list.append(positive_df)
        table_list.append(negative_df)
    print(f"    Added {len(positive_df):,d} positive calls, {len(negative_df):,d} negative calls")

    # EHdn
    if coverage_label != "exome":
        additional_label = ""
        if "downsampled" in results_directory:
            coverage = results_directory.replace("results_for_downsampled_", "").replace("_bam", "")
            additional_label = f"downsampled_to_{coverage}."

        ehdn_table_path = os.path.join(tool_comparison_base_dir, results_directory,
            f"expansion_hunter_denovo/CHM1_CHM13_WGS2.{additional_label}truth_set_EHdn_comparison_table.tsv")

        print(f"Loading {ehdn_table_path}")
        ehdn_df = pd.read_table(ehdn_table_path, dtype=str)
        ehdn_df.loc[:, "Coverage"] = coverage_label
        ehdn_tables.append(ehdn_df)
        print(f"    Added {len(ehdn_df):,d} EHdn calls")


#%%

print(f"Combining {len(alleles_tables)} allele tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.alleles.tsv.gz")
alleles_df = pd.concat(alleles_tables)
alleles_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(alleles_df):,d} rows to {output_path}")
print(alleles_df.groupby(["Coverage", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele: ExpansionHunter - Truth",
    "DiffRepeats: Allele: GangSTR - Truth",
    "DiffRepeats: Allele: HipSTR - Truth",
]])
alleles_df = None

print(f"Combining {len(variants_tables)} variant tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.variants.tsv.gz")
variants_df = pd.concat(variants_tables)
variants_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(variants_df):,d} rows to {output_path}")
print(variants_df.groupby(["Coverage", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele 1: ExpansionHunter - Truth",
    "DiffRepeats: Allele 1: GangSTR - Truth",
    "DiffRepeats: Allele 1: HipSTR - Truth",
    "DiffRepeats: Allele 2: ExpansionHunter - Truth",
    "DiffRepeats: Allele 2: GangSTR - Truth",
    "DiffRepeats: Allele 2: HipSTR - Truth",
]])
variants_df = None

print(f"Combining {len(ehdn_tables)} EHdn tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.EHdn.tsv.gz")
ehdn_df = pd.concat(ehdn_tables)
ehdn_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(ehdn_df):,d} rows to {output_path}")
print(ehdn_df.groupby(["Coverage"]).count()[[
    "LocusId",
]])
ehdn_df = None



#%%