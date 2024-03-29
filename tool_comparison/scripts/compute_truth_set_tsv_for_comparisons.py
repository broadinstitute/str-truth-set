"""This script takes the truth set tsv and the set of negative (non-variant) loci and harmonizes columns across the
2 tables + adds some extra columns that will be useful for tool comparisons.
"""

import argparse
import os
import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

IS_TRUTH_SET_V2_FORMAT = False
HET_or_HOM_or_MULTI_COLUMN = "HET_or_HOM_or_MULTI"
TSV_HEADER = [
    "LocusId", "LocusSize (bp)", "NumRepeatsInReference", "SkippedValidation", "NumRepeatsInT2T",
    "Motif", "CanonicalMotif", "MotifSize",
    "Genotype", "GenotypeConfidenceInterval",
    "NumRepeats: Allele 1", "RepeatSize (bp): Allele 1", "CI start: Allele 1", "CI end: Allele 1", "CI size: Allele 1",
    "NumRepeats: Allele 2", "RepeatSize (bp): Allele 2", "CI start: Allele 2", "CI end: Allele 2", "CI size: Allele 2",
    "DiffFromRefRepeats: Allele 1", "DiffFromRefSize (bp): Allele 1",
    "DiffFromRefRepeats: Allele 2", "DiffFromRefSize (bp): Allele 2",
] + [
    "TruthSetOrNegativeLocus",
    HET_or_HOM_or_MULTI_COLUMN,
    "IsHomRef",
    "IsRef: Allele 1",
    "IsRef: Allele 2",
    "IsPureRepeat",
    "IsFoundInReference",
    "IsMultiallelic",
    "SummaryString",
    "OverlapsIlluminaSTRCatalog: Locus",
    "OverlapsIlluminaSTRCatalog: Motif",
    "OverlapsTRFPureRepeats9bp: Locus",
    "OverlapsTRFPureRepeats9bp: Motif",
    #"OverlapsGangSTRCatalog17: Locus",
    #"OverlapsGangSTRCatalog17: Motif",
    "OverlapsSegDupIntervals",
    "GeneRegionFromGencode_V42",
    "GeneRegionFromMane_V1",
]


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-n", help="Process only the 1st N rows of the truth set TSV. Useful for testing.", type=int)
    p.add_argument("--output-dir", default="./tool_comparison/results/", help="Output directory")
    p.add_argument("truth_set_variants_tsv", help="Path of the truth set variants tsv")
    p.add_argument("negative_loci_tsv", nargs="?", help="Path of negative loci tsv")

    args = p.parse_args()

    if not os.path.isdir(args.output_dir):
        p.error(f"{args.output_dir} doesn't exist")

    for path in args.truth_set_variants_tsv, args.negative_loci_tsv:
        if path and not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def trim_end_column(df):
    """Make sure the locus coordinates span an integer multiple of the motif size. Truncate fractional motif."""
    trimmed_values = df["End1Based"] - (df["End1Based"] - df["Start0Based"]) % df["MotifSize"]
    trimmed_count = sum(df.loc[:, "End1Based"] != trimmed_values)
    if trimmed_count > 0:
        print(f"Trimmed {trimmed_count} out of {len(df)} ({100*trimmed_count/len(df):0.1f}%) locus coordinates")
    df.loc[:, "End1Based"] = trimmed_values


def add_extra_columns(df):
    """Add extra columns to the truth set genotypes so that they have the same columns that will be
    there for tool output genotypes. This will make it easier to compare the two.
    """
    df.loc[:, "LocusSize (bp)"] = df["End1Based"] - df["Start0Based"]

    df.loc[:, "NumRepeats: Allele 1"] = df["NumRepeatsShortAllele"]
    df.loc[:, "RepeatSize (bp): Allele 1"] = df["RepeatSizeShortAllele (bp)"]
    df.loc[:, "CI start: Allele 1"] = df["NumRepeatsShortAllele"]
    df.loc[:, "CI end: Allele 1"] = df["NumRepeatsShortAllele"]
    df.loc[:, "CI size: Allele 1"] = 0

    df.loc[:, "NumRepeats: Allele 2"] = df["NumRepeatsLongAllele"]
    df.loc[:, "RepeatSize (bp): Allele 2"] = df["RepeatSizeLongAllele (bp)"]
    df.loc[:, "CI start: Allele 2"] = df["NumRepeatsLongAllele"]
    df.loc[:, "CI end: Allele 2"] = df["NumRepeatsLongAllele"]
    df.loc[:, "CI size: Allele 2"] = 0

    df.loc[:, "DiffFromRefRepeats: Allele 1"] = df["NumRepeats: Allele 1"] - df["NumRepeatsInReference"]
    df.loc[:, "DiffFromRefSize (bp): Allele 1"] = df["DiffFromRefRepeats: Allele 1"] * df["MotifSize"]
    df.loc[:, "DiffFromRefRepeats: Allele 2"] = df["NumRepeats: Allele 2"] - df["NumRepeatsInReference"]
    df.loc[:, "DiffFromRefSize (bp): Allele 2"] = df["DiffFromRefRepeats: Allele 2"] * df["MotifSize"]

    df.loc[:, "Genotype"] = df["NumRepeats: Allele 1"].astype(str) + "/" + df["NumRepeats: Allele 2"].astype(str)
    df.loc[:, "GenotypeConfidenceInterval"] = df["NumRepeats: Allele 1"].astype(str) + "-" + df["NumRepeats: Allele 1"].astype(str) + "/" + df["NumRepeats: Allele 2"].astype(str) + "-" + df["NumRepeats: Allele 2"].astype(str)

    df.loc[:, "IsRef: Allele 1"] = df["LocusSize (bp)"] == df["RepeatSize (bp): Allele 1"]
    df.loc[:, "IsRef: Allele 2"] = df["LocusSize (bp)"] == df["RepeatSize (bp): Allele 2"]
    df.loc[:, "IsHomRef"] = (
        df["RepeatSize (bp): Allele 1"] == df["RepeatSize (bp): Allele 2"]
    ) & (
        df["LocusSize (bp)"] == df["RepeatSize (bp): Allele 1"]
    )
    error_count = sum(df["IsHomRef"] & (df[HET_or_HOM_or_MULTI_COLUMN] != "HOM"))
    if error_count:
        raise ValueError(f"{error_count} rows have IsHomRef == True, but HET_or_HOM_or_MULTI != HOM")


def main():
    global IS_TRUTH_SET_V2_FORMAT
    global HET_or_HOM_or_MULTI_COLUMN

    args = parse_args()

    # process truth set loci
    truth_set_variants_df = pd.read_table(args.truth_set_variants_tsv, low_memory=False, nrows=args.n, dtype={"Chrom": str})

    IS_TRUTH_SET_V2_FORMAT = "HET_or_HOM_or_MULTI" not in truth_set_variants_df.columns and "HET_or_HOM_or_HEMI_or_MULTI" in truth_set_variants_df.columns
    if IS_TRUTH_SET_V2_FORMAT:
        # update format to match the updated HPRC truth set format
        HET_or_HOM_or_MULTI_COLUMN = "HET_or_HOM_or_HEMI_or_MULTI"
        TSV_HEADER[TSV_HEADER.index("HET_or_HOM_or_MULTI")] = "HET_or_HOM_or_HEMI_or_MULTI"
        # delete "SkippedValidation", "NumRepeatsInT2T" from TSV_HEADER
        del TSV_HEADER[TSV_HEADER.index("SkippedValidation")]
        del TSV_HEADER[TSV_HEADER.index("NumRepeatsInT2T")]

    truth_set_variants_df.loc[:, "TruthSetOrNegativeLocus"] = "TruthSet"
    truth_set_variants_df.loc[:, "Start0Based"] = truth_set_variants_df["Start1Based"] - 1
    truth_set_variants_df.loc[:, "LocusId"] = truth_set_variants_df["LocusId"].str.replace("^chr", "", regex=True)
    trim_end_column(truth_set_variants_df)

    output_path = os.path.basename(args.truth_set_variants_tsv).replace(".tsv", ".for_comparison.tsv")
    write_to_tsv(truth_set_variants_df, os.path.join(args.output_dir, output_path))

    # process negative loci
    if args.negative_loci_tsv:
        negative_loci_df = pd.read_table(args.negative_loci_tsv, low_memory=False, nrows=args.n, dtype={"Chrom": str})
        trim_end_column(negative_loci_df)
        negative_loci_df.loc[:, "TruthSetOrNegativeLocus"] = "NegativeLocus"
        negative_loci_df.loc[:, "MotifSize"] = negative_loci_df["Motif"].str.len()
        negative_loci_df.loc[:, "Start1Based"] = negative_loci_df["Start0Based"] + 1
        negative_loci_df.loc[:, "Locus"] = negative_loci_df["Chrom"].astype(str) + ":" + \
                                           negative_loci_df["Start1Based"].astype(str) + "-" + \
                                           negative_loci_df["End1Based"].astype(str)
        negative_loci_df.loc[:, "Chrom"] = negative_loci_df["Chrom"].str.replace("^chr", "", regex=True)
        negative_loci_df.loc[:, "LocusId"] = negative_loci_df[["Chrom", "Start0Based", "End1Based", "Motif"]].astype(str).apply("-".join, axis=1)
        negative_loci_df.loc[:, "RepeatSizeShortAllele (bp)"] = negative_loci_df["End1Based"] - negative_loci_df["Start0Based"]
        negative_loci_df.loc[:, "RepeatSizeLongAllele (bp)"] = negative_loci_df["RepeatSizeShortAllele (bp)"]
        negative_loci_df.loc[:, "NumRepeatsInReference"] = (
                (negative_loci_df["End1Based"] - negative_loci_df["Start0Based"])/negative_loci_df["MotifSize"]
        ).astype(int)
        if not IS_TRUTH_SET_V2_FORMAT:
            negative_loci_df.loc[:, "NumRepeatsInT2T"] = ""
            negative_loci_df.loc[:, "SkippedValidation"] = True
        negative_loci_df.loc[:, "NumRepeatsShortAllele"] = negative_loci_df.loc[:, "NumRepeatsInReference"]
        negative_loci_df.loc[:, "NumRepeatsLongAllele"] = negative_loci_df.loc[:, "NumRepeatsInReference"]
        negative_loci_df.loc[:, HET_or_HOM_or_MULTI_COLUMN] = "HOM"
        negative_loci_df.loc[:, "IsPureRepeat"] = True
        negative_loci_df.loc[:, "IsFoundInReference"] = True
        negative_loci_df.loc[:, "IsMultiallelic"] = False
        negative_loci_df.loc[:, "SummaryString"] = ""
        negative_loci_df.loc[:, "CanonicalMotif"] = negative_loci_df.Motif.apply(
            lambda m: compute_canonical_motif(m, include_reverse_complement=True))
        negative_loci_df.loc[:, "SummaryString"] = "RU" + negative_loci_df["MotifSize"].astype(str) + ":" + negative_loci_df["Motif"] + ":NEGATIVE-HOM-REF:" + negative_loci_df["NumRepeatsInReference"].astype(str)

        output_path = os.path.basename(args.negative_loci_tsv).replace(".tsv", ".for_comparison.tsv")
        write_to_tsv(negative_loci_df, os.path.join(args.output_dir, output_path))


def write_to_tsv(df, output_path):
    add_extra_columns(df)
    df = df[TSV_HEADER]
    print(f"Output columns for {output_path}:")
    for column in sorted(df.columns):
        is_nan_count = sum(pd.isna(df[column]))
        print(f"\t{is_nan_count:7,d} of {len(df):,d} ({100*is_nan_count/len(df):5.1f}%) NaN values in {column}")

    df.to_csv(output_path, header=True, index=False, sep="\t")

    print(f"Wrote {len(df):,d} rows to {output_path}")


if __name__ == "__main__":
    main()
