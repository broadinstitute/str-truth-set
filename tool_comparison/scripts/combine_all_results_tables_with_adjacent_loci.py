import os
import pandas as pd

tool_comparison_base_dir = "tool_comparison"


def locus_id_to_region(locus_id):
    locus_id = locus_id.split(";")[0]
    chrom, start_0based, end, motif = locus_id.split("-")
    return f"{chrom}:{int(start_0based) + 1}-{end}"


variants_tables = []
alleles_tables = []
for adjacent_loci_label, results_directory in [
    ("no_adjacent_loci", "results"),
    ("max_dist_10bp", "results_with_adjacent_loci__max_dist_10bp"),
    ("max_dist_24bp", "results_with_adjacent_loci__max_dist_24bp"),
    ("max_dist_50bp", "results_with_adjacent_loci__max_dist_50bp"),
]:
    # ExpansionHunter
    print(f"Processing {results_directory}")
    for table_list, filename_suffix in [
        (variants_tables, ".variants"),
        (alleles_tables, ".alleles"),
    ]:
        positive_df_path = os.path.join(
            tool_comparison_base_dir, results_directory, f"STR_truth_set.v1.for_comparison{filename_suffix}.tsv")
        print(f"Loading {positive_df_path}")
        positive_df = pd.read_table(positive_df_path, dtype=str)

        positive_df.loc[:, "Locus"] = positive_df["LocusId"].apply(locus_id_to_region)
        positive_df.loc[:, "Coverage"] = "40x"
        positive_df.loc[:, "AdjacentLociLabel"] = adjacent_loci_label
        positive_df.loc[:, "PositiveOrNegative"] = "positive"
        positive_df = positive_df[
            [c for c in positive_df.columns if not any(k in c.lower() for k in ("gangstr", "hipstr"))]
        ]
        table_list.append(positive_df)
    print(f"    Added {len(positive_df):,d} positive calls")


#%%

print(f"Combining {len(alleles_tables)} allele tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.with_adjacent_loci.alleles.tsv.gz")
alleles_df = pd.concat(alleles_tables)
alleles_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(alleles_df):,d} rows to {output_path}")
print(alleles_df.groupby(["AdjacentLociLabel", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele: ExpansionHunter - Truth",
]])
alleles_df = None

print(f"Combining {len(variants_tables)} variant tables..")
output_path = os.path.join(tool_comparison_base_dir, "combined.results.with_adjacent_loci.variants.tsv.gz")
variants_df = pd.concat(variants_tables)
variants_df.to_csv(output_path, sep="\t", header=True, index=False)
print(f"Wrote {len(variants_df):,d} rows to {output_path}")
print(variants_df.groupby(["AdjacentLociLabel", "PositiveOrNegative"]).count()[[
    "LocusId",
    "DiffRepeats: Allele 1: ExpansionHunter - Truth",
    "DiffRepeats: Allele 2: ExpansionHunter - Truth",
]])
variants_df = None



#%%