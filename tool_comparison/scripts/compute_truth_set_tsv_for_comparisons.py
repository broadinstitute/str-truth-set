import argparse
import os
import pandas as pd

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

TSV_HEADER = [
    "LocusId", "LocusSize (bp)", "Motif", "MotifSize",
    "Genotype", "GenotypeConfidenceInterval",
    "NumRepeats: Allele 1", "RepeatSize (bp): Allele 1", "CI start: Allele 1", "CI end: Allele 1", "CI size: Allele 1",
    "NumRepeats: Allele 2", "RepeatSize (bp): Allele 2", "CI start: Allele 2", "CI end: Allele 2", "CI size: Allele 2",
    "DiffFromRefRepeats: Allele 1", "DiffFromRefSize (bp): Allele 1",
    "DiffFromRefRepeats: Allele 2", "DiffFromRefSize (bp): Allele 2",
] + [
    "HET_or_HOM",
    "IsHomRef",
    "IsFoundInReference",
    "SummaryString",
    "IsMultiallelic",
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
    p.add_argument("--output-dir", default="./tool_comparison/results/", help="Output directory")
    p.add_argument("truth_set_variants_tsv", help="Path of the truth set variants tsv")
    p.add_argument("negative_loci_tsv", help="Path of negative loci tsv")

    args = p.parse_args()

    if not os.path.isdir(args.output_dir):
        p.error(f"{args.output_dir} doesn't exist")

    for path in args.truth_set_variants_tsv, args.negative_loci_tsv:
        if path and not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def trim_end_column(df):
    trimmed_values = df["End1Based"] - (df["End1Based"] - df["Start0Based"]) % df["MotifSize"]
    trimmed_count = sum(df.loc[:, "End1Based"] != trimmed_values)
    if trimmed_count > 0:
        print(f"Trimmed {trimmed_count} out of {len(df)} ({100*trimmed_count/len(df):0.1f}%) locus coordinates")
    df.loc[:, "End1Based"] = trimmed_values


def add_extra_columns(df):
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
    df.loc[:, "IsHomRef"] = (
        df["RepeatSize (bp): Allele 1"] == df["RepeatSize (bp): Allele 2"]
    ) & (
        df["LocusSize (bp)"] == df["RepeatSize (bp): Allele 1"]
    )
    error_count = sum(df["IsHomRef"] & (df["HET_or_HOM"] != "HOM"))
    if error_count:
        raise ValueError(f"{error_count} rows have IsHomRef == True, but HET_or_HOM != HOM")


def main():
    args = parse_args()

    truth_set_variants_df = pd.read_table(args.truth_set_variants_tsv, low_memory=False)
    truth_set_variants_df.loc[:, "Start0Based"] = truth_set_variants_df["Start1Based"] - 1
    truth_set_variants_df.loc[:, "TruthSetOrNegativeLocus"] = "TruthSet"
    truth_set_variants_df.loc[:, "LocusId"] = truth_set_variants_df["LocusId"].str.replace("^chr", "", regex=True)

    output_path = os.path.basename(args.truth_set_variants_tsv).replace(".tsv", ".for_comparison.tsv")
    write_to_tsv(truth_set_variants_df, os.path.join(args.output_dir, output_path))

    negative_loci_df = pd.read_table(args.negative_loci_tsv, low_memory=False)
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
    negative_loci_df.loc[:, "NumRepeatsShortAllele"] = negative_loci_df.loc[:, "NumRepeatsInReference"]
    negative_loci_df.loc[:, "NumRepeatsLongAllele"] = negative_loci_df.loc[:, "NumRepeatsInReference"]
    negative_loci_df.loc[:, "HET_or_HOM"] = "HOM"
    negative_loci_df.loc[:, "IsFoundInReference"] = True
    negative_loci_df.loc[:, "SummaryString"] = ""
    negative_loci_df.loc[:, "IsMultiallelic"] = False
    negative_loci_df.loc[:, "CanonicalMotif"] = negative_loci_df.Motif.apply(
        lambda m: compute_canonical_motif(m, include_reverse_complement=True))
    negative_loci_df.loc[:, "SummaryString"] = "RU" + negative_loci_df["MotifSize"].astype(str) + ":" + negative_loci_df["Motif"] + ":NEGATIVE-HOM-REF:" + negative_loci_df["NumRepeatsInReference"].astype(str)

    output_path = os.path.basename(args.negative_loci_tsv).replace(".tsv", ".for_comparison.tsv")
    write_to_tsv(negative_loci_df, os.path.join(args.output_dir, output_path))


def write_to_tsv(df, output_path):
    trim_end_column(df)
    add_extra_columns(df)

    df = df[TSV_HEADER]
    df.to_csv(output_path, header=True, index=False, sep="\t")
    print("Output columns:")
    for column in sorted(df.columns):
        is_nan_count = sum(pd.isna(df[column]))
        print(f"\t{is_nan_count:7,d} of {len(df):,d} ({100*is_nan_count/len(df):5.1f}%) NaN values in {column}")

    print(f"Wrote {len(df):,d} rows to {output_path}")


if __name__ == "__main__":
    main()
