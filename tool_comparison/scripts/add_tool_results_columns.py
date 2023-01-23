import argparse
import os
import pandas as pd
from pprint import pprint
import re

TOOL_DF_COLUMNS_TO_KEEP = [
    "LocusId",
    "Coverage",
    "Motif",
    "MotifSize",
    "Genotype",
    "GenotypeConfidenceInterval",

    "IsHomRef",
    "IsRef: Allele 1",
    "IsRef: Allele 2",

    "NumRepeats: Allele 1",
    "RepeatSize (bp): Allele 1",
    "CI start: Allele 1",
    "CI end: Allele 1",
    "CI size: Allele 1",
    "NumRepeats: Allele 2",
    "RepeatSize (bp): Allele 2",
    "CI start: Allele 2",
    "CI end: Allele 2",
    "CI size: Allele 2",

    "DiffFromRefRepeats: Allele 1",
    "DiffFromRefRepeats: Allele 2",
    "DiffFromRefSize (bp): Allele 1",
    "DiffFromRefSize (bp): Allele 2",
]

EH_AND_GANGSTR_COLUMNS = [
    "NumReadsTotal",
    "NumSpanningReads",
    "NumFlankingReads",
    "NumInrepeatReads",
    "FractionOfReadsThatSupportsGenotype: Allele 1",
    "FractionOfReadsThatSupportsGenotype: Allele 2",
    "NumReadsTotalThatSupportGenotype: Allele 1",
    "NumReadsTotalThatSupportGenotype: Allele 2",
    "NumSpanningReadsThatSupportGenotype: Allele 1",
    "NumSpanningReadsThatSupportGenotype: Allele 2",
]

MERGE_KEY_COLUMNS = ["LocusId", "Motif", "MotifSize"]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--verbose", action="store_true", help="Whether to print additional info about input and output columns.")
    p.add_argument("--output-tsv", help="Output path of combined tsv file")
    p.add_argument("--tool", choices={"ExpansionHunter", "GangSTR", "HipSTR"})
    p.add_argument("tool_results_tsv", help="Path of the tool results combined tsv file.")
    p.add_argument("truth_set_or_negative_loci_tsv", help="Path of the truth set or negative_loci tsv")

    args = p.parse_args()

    for path in args.tool_results_tsv, args.truth_set_or_negative_loci_tsv:
        if path and not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def main():
    args = parse_args()

    tool_df_columns_to_keep = list(TOOL_DF_COLUMNS_TO_KEEP)
    if args.tool == "ExpansionHunter":
        tool_df_columns_to_keep += EH_AND_GANGSTR_COLUMNS
        tool_df_columns_to_keep += ["NumAllelesSupportedTotal"]
    elif args.tool == "GangSTR":
        tool_df_columns_to_keep += EH_AND_GANGSTR_COLUMNS
        tool_df_columns_to_keep += ["Q"]
    elif args.tool == "HipSTR":
        tool_df_columns_to_keep += ["Q", "DP", "AB", "FS", "DFLANKINDEL", "DSTUTTER"]
    else:
        raise ValueError(f"Unexpected tool: {args.tool}")
    
    truth_set_df = pd.read_table(args.truth_set_or_negative_loci_tsv)

    tool_df = pd.read_table(args.tool_results_tsv)
    tool_df.rename(columns={
        "RepeatUnit": "Motif",
        "RepeatUnitLength": "MotifSize",
        "Repeat Size (bp): Allele 1": "RepeatSize (bp): Allele 1",
        "Repeat Size (bp): Allele 2": "RepeatSize (bp): Allele 2",
        "Num Repeats: Allele 1": "NumRepeats: Allele 1",
        "Num Repeats: Allele 2": "NumRepeats: Allele 2",
    }, inplace=True)

    tool_df.loc[:, "ReferenceRegion"] = tool_df["ReferenceRegion"].str.replace("^chr", "", regex=True)
    tool_df.loc[:, "LocusId"] = tool_df["LocusId"].str.replace("^chr", "", regex=True)

    def split_reference_region(row):
        result = re.split("[:-]", row["ReferenceRegion"])
        if len(result) != 3:
            raise ValueError(f"Unexpected ReferenceRegion: " + str(row["ReferenceRegion"]))
        return result

    tool_df[["Chrom", "Start0Based", "End1Based"]] = tool_df.apply(split_reference_region, axis=1, result_type="expand")
    tool_df.loc[:, "Start0Based"] = tool_df["Start0Based"].astype(int)
    tool_df.loc[:, "End1Based"] = tool_df["End1Based"].astype(int)

    tool_df.loc[:, "LocusSize (bp)"] = tool_df["End1Based"] - tool_df["Start0Based"]
    tool_df.loc[:, "NumRepeatsInReference"] = (tool_df["LocusSize (bp)"]/tool_df["MotifSize"]).astype(int)
    tool_df.loc[:, "DiffFromRefRepeats: Allele 1"] = tool_df["NumRepeats: Allele 1"] - tool_df["NumRepeatsInReference"]
    tool_df.loc[:, "DiffFromRefRepeats: Allele 2"] = tool_df["NumRepeats: Allele 2"] - tool_df["NumRepeatsInReference"]
    tool_df.loc[:, "DiffFromRefSize (bp): Allele 1"] = tool_df["DiffFromRefRepeats: Allele 1"] * tool_df["MotifSize"]
    tool_df.loc[:, "DiffFromRefSize (bp): Allele 2"] = tool_df["DiffFromRefRepeats: Allele 2"] * tool_df["MotifSize"]

    tool_df.loc[:, "IsRef: Allele 1"] = tool_df["LocusSize (bp)"].astype(int) == tool_df["RepeatSize (bp): Allele 1"].astype(int)
    tool_df.loc[:, "IsRef: Allele 2"] = tool_df["LocusSize (bp)"].astype(int) == tool_df["RepeatSize (bp): Allele 2"].astype(int)
    tool_df.loc[:, "IsHomRef"] = (
         tool_df["RepeatSize (bp): Allele 1"] == tool_df["RepeatSize (bp): Allele 2"]
    ) & (
        tool_df["LocusSize (bp)"] == tool_df["RepeatSize (bp): Allele 1"]
    )

    if args.verbose:
        print("="*100)
        for column in sorted(tool_df.columns):
            is_nan_count = sum(pd.isna(tool_df[column]))
            print(f"\t{is_nan_count:7,d} of {len(tool_df):,d} ({100*is_nan_count/len(tool_df):5.1f}%) NaN values in {column}")

        print("="*100)
        print("set(truth_set_df.columns) & set(tool_df.columns)")
        pprint(sorted(set(truth_set_df.columns) & set(tool_df.columns)))
        print("="*100)
        print("set(truth_set_df.columns) - set(tool_df.columns)")
        pprint(sorted(set(truth_set_df.columns) - set(tool_df.columns)))
        print("="*100)
        print("set(tool_df.columns) - set(truth_set_df.columns)")
        pprint(sorted(set(tool_df.columns) - set(truth_set_df.columns)))
        print("="*100)
        print("set(tool_df.columns) - set(columns to keep)")
        pprint(sorted(set(tool_df.columns) - set(tool_df_columns_to_keep)))

    tool_df = tool_df[tool_df_columns_to_keep]
    for tool_column in set(tool_df.columns) - set(truth_set_df.columns):
        tool_df.rename(columns={tool_column: f"{tool_column}: {args.tool}"}, inplace=True)

    df_merged = pd.merge(
        truth_set_df,
        tool_df,
        how="left",
        on=MERGE_KEY_COLUMNS,
        suffixes=(f": Truth", f": {args.tool}"))

    if len(df_merged) > len(truth_set_df):
        print(f"WARNING: merged column has {len(df_merged) - len(truth_set_df):,d} more rows than truth set.")
    elif len(df_merged) == len(truth_set_df):
        print(f"Merged table has the same number of rows ({len(truth_set_df):,d}) as the truth set.")
    else:
        raise ValueError("Merged table has fewer rows than the results table.")

    if not args.output_tsv:
        args.output_tsv = args.truth_set_or_negative_loci_tsv.replace(".tsv", f".with_{args.tool}_results.tsv")

    df_merged.to_csv(args.output_tsv, index=False, header=True, sep="\t")

    print("Output columns:")
    for column in sorted(df_merged.columns):
        is_nan_count = sum(pd.isna(df_merged[column]))
        print(f"\t{is_nan_count:7,d} of {len(df_merged):,d} ({100*is_nan_count/len(df_merged):5.1f}%) NaN values in {column}")
    print(f"Wrote {len(df_merged):,d} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()
