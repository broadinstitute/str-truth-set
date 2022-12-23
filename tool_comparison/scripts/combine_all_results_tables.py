import glob
import os
import pandas as pd

tool_comparison_base_dir = "tool_comparison"

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
    print(f"Processing {results_directory}")
    for table_list, filename_suffix in [
        (variants_tables, ".variants"),
        (alleles_tables, ".alleles"),
    ]:
        positive_df = pd.read_table(os.path.join(
            tool_comparison_base_dir, results_directory, f"STR_truthset.v1.for_comparison{filename_suffix}.tsv"),
            dtype=str,
        )
        negative_df = pd.read_table(os.path.join(
            tool_comparison_base_dir, results_directory, f"negative_loci.for_comparison{filename_suffix}.tsv"),
            dtype=str,
        )

        positive_df.loc[:, "coverage"] = coverage_label
        negative_df.loc[:, "coverage"] = coverage_label
        positive_df.loc[:, "PositiveOrNegative"] = "positive"
        negative_df.loc[:, "PositiveOrNegative"] = "negative"

        table_list.append(positive_df)
        table_list.append(negative_df)
    print(f"    Added {len(positive_df):,d} positive calls, {len(negative_df):,d} negative calls")

    if coverage_label != "exome":
        ehdn_table_wildcard_path = os.path.join(tool_comparison_base_dir, results_directory,
            f"expansion_hunter_denovo/CHM1_CHM13_WGS2.*truth_set_EHdn_comparison_table.tsv")
        ehdn_table_paths = glob.glob(ehdn_table_wildcard_path)

        if len(ehdn_table_paths) != 1:
            raise ValueError(f"Found {len(ehdn_table_paths)} EHdn tables in {ehdn_table_wildcard_path}")

        ehdn_df = pd.read_table(ehdn_table_paths[0], dtype=str)
        ehdn_df.loc[:, "coverage"] = coverage_label
        ehdn_tables.append(ehdn_df)
        print(f"    Added {len(ehdn_df):,d} EHdn calls")


#%%

print(f"Combining {len(alleles_tables)} allele tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.alleles.tsv")
alleles_df = pd.concat(alleles_tables)
alleles_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(alleles_df):,d} rows to {output_path}")
print(alleles_df.groupby(["coverage", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele: ExpansionHunter - Truth",
    "DiffRepeats: Allele: GangSTR - Truth",
]])
alleles_df = None

print(f"Combining {len(variants_tables)} variant tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.variants.tsv")
variants_df = pd.concat(variants_tables)
variants_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(variants_df):,d} rows to {output_path}")
print(variants_df.groupby(["coverage", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele 1: ExpansionHunter - Truth",
    "DiffRepeats: Allele 1: GangSTR - Truth",
    "DiffRepeats: Allele 2: ExpansionHunter - Truth",
    "DiffRepeats: Allele 2: GangSTR - Truth",
]])
variants_df = None

print(f"Combining {len(ehdn_tables)} EHdn tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.EHdn.tsv")
edhn_df = pd.concat(ehdn_tables)
edhn_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(edhn_df):,d} rows to {output_path}")
print(edhn_df.groupby(["coverage"]).count()[[
    "LocusId",
]])
edhn_df = None



#%%