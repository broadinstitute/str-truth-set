"""Add sequence-accuracy columns to a tool-vs-truth comparison table.

The accuracy plots score a tool by how close its repeat *count* is to the truth. This script scores the tools that
report an actual allele *sequence* (TRGT, ATaRVa, HipSTR) by the edit distance between the sequence they report and
the assembly-derived truth allele sequence, which also captures the substitutions and interruptions that a repeat
count can't express.

It joins two side tables onto the ...with_{tool}_results.tsv.gz table produced by add_tool_results_columns.py:

    --truth-set-genotypes   {sample}.tandem_repeat_genotypes.tsv.gz, whose Allele1Sequence / Allele2Sequence columns
                            hold the assembly's two haplotype sequences spanning the locus interval
    --allele-sequences      {prefix}.allele_sequences.tsv.gz from str_analysis.extract_allele_sequences_from_vcf,
                            which holds the tool's two allele sequences over the same interval

and adds 4 numeric columns:

    SequenceEditDistance: Allele 1: {tool}             Levenshtein distance in bp
    SequenceEditDistance: Allele 2: {tool}
    SequenceEditDistanceNormalized: Allele 1: {tool}   the same distance / max(1, len(truth allele sequence))
    SequenceEditDistanceNormalized: Allele 2: {tool}

The "Allele N: " naming is required: add_concordance_columns.write_alleles_table() raises a ValueError on any column
that contains "Allele 1" without the exact substring "Allele 1: ", and it melts these columns into
"SequenceEditDistance: Allele: {tool}" for free.

Every sequence column is dropped before the table is written, so the sequences (which run to kilobases for VNTRs)
never reach the combined variants/alleles tables.

An allele scores NA rather than a distance whenever either side's sequence is unavailable, so an unavailable sequence
is never silently scored as a perfect match against an empty string. See pair_truth_sequences() and
read_tool_allele_sequences() for the two presence rules.
"""

import argparse
import os
import pandas as pd

from rapidfuzz.distance import Levenshtein

from str_analysis.extract_allele_sequences_from_vcf import CALLED, EXTRACTION_OK

# Tools scored by this benchmark. Must match run_tools/run_genotyping_tools.py's SEQUENCE_ACCURACY_TOOLS --
# TRGTv3 is scoped out of v1 there, so it's excluded here too even though its VCF format is TRGTv5-compatible.
SEQUENCE_ACCURACY_TOOLS = ["TRGTv5", "ATaRVa", "HipSTR"]

TRUTH_SEQUENCE_COLUMNS = ["TruthAlleleSequence: Allele 1", "TruthAlleleSequence: Allele 2"]
TOOL_SEQUENCE_COLUMNS = ["ToolAlleleSequence: Allele 1", "ToolAlleleSequence: Allele 2"]


def sequence_or_none(sequence):
    """Return the sequence, or None if it isn't a string (a left-join miss arrives as NaN, not None)."""
    return sequence if isinstance(sequence, str) else None


def pair_truth_sequences(allele_1_sequence, allele_2_sequence, short_allele_size_bp, long_allele_size_bp):
    """Return the truth set's two allele sequences ordered short-first, or (None, None) if they aren't usable.

    filter_vcf_to_tandem_repeats writes Allele1Sequence / Allele2Sequence per haplotype (haplotype 0 and haplotype 1)
    rather than sorted by size, while every downstream table uses Allele 1 = short and Allele 2 = long. Sorting the
    two sequences by length reproduces that ordering: across the whole HG002 truth set (583,174 loci) the sorted
    lengths equal (RepeatSizeShortAlleleBp, RepeatSizeLongAlleleBp) for every locus where both sequences are stored.

    That length check doubles as the presence rule, which matters because an empty sequence means two different
    things. It is a real zero-length allele exactly when its reported size is 0, and is an unavailable sequence
    otherwise -- so an unavailable sequence is never scored as if the tool had to reproduce "".

    The one systematic case where a sequence is missing is HEMI loci (male chrX/chrY outside the PAR): the truth set
    stores one real sequence plus an empty string, but duplicates the single allele size into both RepeatSize columns
    (19,050 of HG002's 19,098 HEMI loci). Those are read as a duplicated diploid genotype, matching what the truth set
    already does for the sizes and what add_tool_results_columns.py does for a tool's haploid calls. Without this,
    HEMI loci would contribute one scored allele against the tool's two.

    Args:
        allele_1_sequence (str): the haplotype 0 sequence ("" when the truth set stored none).
        allele_2_sequence (str): the haplotype 1 sequence ("" when the truth set stored none).
        short_allele_size_bp (float): the locus's RepeatSizeShortAlleleBp.
        long_allele_size_bp (float): the locus's RepeatSizeLongAlleleBp.

    Returns:
        tuple: (short allele sequence, long allele sequence), each a string or None.
    """
    short_sequence, long_sequence = sorted([allele_1_sequence, allele_2_sequence], key=len)
    short_matches = len(short_sequence) == short_allele_size_bp
    long_matches = len(long_sequence) == long_allele_size_bp
    if short_matches and long_matches:
        return short_sequence, long_sequence
    if short_allele_size_bp == long_allele_size_bp:
        # HEMI: a single allele size duplicated across both columns, with only one sequence stored
        if long_matches:
            return long_sequence, long_sequence
        if short_matches:
            return short_sequence, short_sequence
    return None, None


def compute_sequence_edit_distances(truth_allele_1, truth_allele_2, tool_allele_1, tool_allele_2):
    """Pair the truth and tool alleles and return the edit distance for each.

    Both sides arrive size-sorted (Allele 1 is the shorter allele), so the index-wise pairing is the size-sorted one.
    When either side's two alleles happen to have the same length, both assignments are equally size-sorted, so both
    are evaluated and the one with the smaller total distance wins; an exact tie keeps the index-wise pairing, which
    makes the result deterministic.

    Args:
        truth_allele_1 (str): the shorter truth allele sequence, or None when it isn't available.
        truth_allele_2 (str): the longer truth allele sequence, or None.
        tool_allele_1 (str): the shorter tool allele sequence, or None when the tool didn't call it.
        tool_allele_2 (str): the longer tool allele sequence, or None.

    Returns:
        tuple: (distance for Allele 1, distance for Allele 2), each an int or None when that allele can't be scored.
    """
    truth_allele_1 = sequence_or_none(truth_allele_1)
    truth_allele_2 = sequence_or_none(truth_allele_2)
    tool_allele_1 = sequence_or_none(tool_allele_1)
    tool_allele_2 = sequence_or_none(tool_allele_2)

    distances = (
        None if truth_allele_1 is None or tool_allele_1 is None else Levenshtein.distance(truth_allele_1, tool_allele_1),
        None if truth_allele_2 is None or tool_allele_2 is None else Levenshtein.distance(truth_allele_2, tool_allele_2),
    )

    if None in (truth_allele_1, truth_allele_2, tool_allele_1, tool_allele_2):
        return distances
    if len(truth_allele_1) != len(truth_allele_2) and len(tool_allele_1) != len(tool_allele_2):
        return distances

    swapped = (
        Levenshtein.distance(truth_allele_1, tool_allele_2),
        Levenshtein.distance(truth_allele_2, tool_allele_1),
    )
    return swapped if sum(swapped) < sum(distances) else distances


def read_truth_allele_sequences(truth_set_genotypes_tsv):
    """Read the truth set's per-allele sequences, ordered short-first by pair_truth_sequences().

    Args:
        truth_set_genotypes_tsv (str): path of {sample}.tandem_repeat_genotypes.tsv.gz.

    Returns:
        pd.DataFrame: LocusId (with any leading "chr" stripped, matching the comparison tables) plus the two
        TRUTH_SEQUENCE_COLUMNS.
    """
    df = pd.read_table(
        truth_set_genotypes_tsv,
        usecols=["LocusId", "Allele1Sequence", "Allele2Sequence",
                 "RepeatSizeShortAlleleBp", "RepeatSizeLongAlleleBp"],
        dtype={"LocusId": str, "Allele1Sequence": str, "Allele2Sequence": str})
    print(f"Read {len(df):,d} rows from {truth_set_genotypes_tsv}")

    sequences = pd.DataFrame([
        pair_truth_sequences(allele_1_sequence, allele_2_sequence, short_allele_size_bp, long_allele_size_bp)
        for allele_1_sequence, allele_2_sequence, short_allele_size_bp, long_allele_size_bp in zip(
            # ANALYSIS_OK[imputation]: "" is a placeholder for len()/sorting only, never a claim that the allele
            # is zero-length. pair_truth_sequences() decides availability from the reported RepeatSize*AlleleBp
            # values rather than from whether the sequence field was empty, so a missing sequence still returns
            # None (scores NA) instead of being scored as a zero-length allele.
            df["Allele1Sequence"].fillna(""), df["Allele2Sequence"].fillna(""),
            df["RepeatSizeShortAlleleBp"], df["RepeatSizeLongAlleleBp"])
    ], index=df.index, columns=TRUTH_SEQUENCE_COLUMNS)
    sequences.loc[:, "LocusId"] = df["LocusId"].str.replace("^chr", "", regex=True)

    unavailable = int(sequences[TRUTH_SEQUENCE_COLUMNS[0]].isna().sum())
    print(f"{unavailable:,d} of {len(sequences):,d} truth loci have no usable allele sequence and will score NA")

    return sequences[["LocusId"] + TRUTH_SEQUENCE_COLUMNS]


def read_tool_allele_sequences(allele_sequences_tsv):
    """Read the extractor's table, keeping only the allele sequences that can be scored.

    An allele is scored only when the locus extracted cleanly (ExtractionStatus == "ok") and that allele was actually
    called (AlleleStatus == "called"). Under those two conditions an empty sequence is a genuine zero-length allele
    (ATaRVa writes <DEL> for these) and gets a real distance; every other empty value becomes None.

    Args:
        allele_sequences_tsv (str): path of {prefix}.allele_sequences.tsv.gz.

    Returns:
        pd.DataFrame: LocusId, the two TOOL_SEQUENCE_COLUMNS, and ExtractionStatus / AlleleStatus columns, which the
        caller uses for the QC printout and then drops.
    """
    # keep_default_na=False so that an empty AlleleSequence stays "" (a zero-length allele) instead of becoming NaN
    df = pd.read_table(allele_sequences_tsv, dtype=str, keep_default_na=False)
    print(f"Read {len(df):,d} rows from {allele_sequences_tsv}")

    for allele in "Allele 1", "Allele 2":
        df[f"ToolAlleleSequence: {allele}"] = [
            sequence if extraction_status == EXTRACTION_OK and allele_status == CALLED else None
            for sequence, allele_status, extraction_status in zip(
                df[f"AlleleSequence: {allele}"], df[f"AlleleStatus: {allele}"], df["ExtractionStatus"])
        ]

    return df[["LocusId", "ExtractionStatus", "AlleleStatus: Allele 1", "AlleleStatus: Allele 2"]
              + TOOL_SEQUENCE_COLUMNS]


def print_qc(df, tool, tool_sequences_df, rows_in):
    """Print the QC block: how the two joins landed, why alleles went to NA, and the distance distribution."""
    print("=" * 100)
    print(f"Sequence accuracy QC for {tool}:")
    print(f"{rows_in:12,d} rows in the comparison table")

    joined = int(df["ExtractionStatus"].notna().sum())
    print(f"{joined:12,d} rows joined to an {tool} allele-sequence row "
          f"({rows_in - joined:,d} loci {tool} produced no record for)")
    unjoined_tool_loci = len(set(tool_sequences_df["LocusId"]) - set(df["LocusId"]))
    print(f"{unjoined_tool_loci:12,d} {tool} allele-sequence rows have no matching row in the comparison table")

    for status, count in sorted(df["ExtractionStatus"].value_counts().items()):
        print(f"{count:12,d} loci with ExtractionStatus == {status}")
    for allele in "Allele 1", "Allele 2":
        for status, count in sorted(df[f"AlleleStatus: {allele}"].value_counts().items()):
            print(f"{count:12,d} {allele} with AlleleStatus == {status}")

    print(f"{int(df['TruthAlleleSequence: Allele 1'].isna().sum()):12,d} rows have no usable truth allele sequence")

    for allele in "Allele 1", "Allele 2":
        # ANALYSIS_OK[imputation]: this feeds a printed QC summary only, never the output table. The NA alleles are
        # already counted by the status lines above, and including them here would dilute the "% exactly right"
        # figure with alleles that were never scored.
        distances = pd.to_numeric(df[f"SequenceEditDistance: {allele}: {tool}"], errors="coerce").dropna()
        if len(distances) == 0:
            print(f"             No {allele} distances were computed")
            continue
        print(f"{len(distances):12,d} {allele} distances: {100 * (distances == 0).mean():.1f}% exactly right, "
              f"median {distances.median():.0f}bp, mean {distances.mean():.1f}bp, "
              f"90th percentile {distances.quantile(0.9):.0f}bp, max {distances.max():.0f}bp")
    print("=" * 100)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tool", choices=SEQUENCE_ACCURACY_TOOLS, required=True,
                   help="Which tool's results are in the input table.")
    p.add_argument("--truth-set-genotypes", required=True,
                   help="Path of {sample}.tandem_repeat_genotypes.tsv.gz (the truth set, which carries the "
                        "assembly's per-haplotype allele sequences).")
    p.add_argument("--allele-sequences", required=True,
                   help="Path of {prefix}.allele_sequences.tsv.gz from str_analysis.extract_allele_sequences_from_vcf.")
    p.add_argument("--output-tsv", help="Output path. Defaults to overwriting the input table in place, since the "
                                        "only change is the 4 added columns.")
    p.add_argument("combined_tsv", help="Path of the ...with_{tool}_results.tsv.gz table from "
                                        "add_tool_results_columns.py.")

    args = p.parse_args()

    for path in args.combined_tsv, args.truth_set_genotypes, args.allele_sequences:
        if not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def main():
    args = parse_args()

    # dtype=str / keep_default_na=False so the columns this script doesn't touch round-trip through the file
    # unchanged (no float reformatting, no "NA"-like strings turning into NaN)
    df = pd.read_table(args.combined_tsv, dtype=str, keep_default_na=False)
    rows_in = len(df)
    print(f"Read {rows_in:,d} rows from {args.combined_tsv}")

    truth_sequences_df = read_truth_allele_sequences(args.truth_set_genotypes)
    tool_sequences_df = read_tool_allele_sequences(args.allele_sequences)

    for side_df in truth_sequences_df, tool_sequences_df:
        # A handful of loci appear twice in the truth set (3 of HG002's 583,174, as byte-identical duplicate rows),
        # and every catalog and tool VCF derived from it inherits them. Both copies carry the same sequences, so
        # keeping one is lossless - but leaving them in would multiply those rows in the merge below.
        duplicate_count = int(side_df["LocusId"].duplicated().sum())
        if duplicate_count:
            print(f"Dropping {duplicate_count:,d} duplicate LocusId rows before merging")
            side_df = side_df.drop_duplicates(subset="LocusId", keep="first")
        # validate="many_to_one" makes any remaining row multiplication an error rather than a silent miscount
        df = df.merge(side_df, how="left", on="LocusId", validate="many_to_one")
        assert len(df) == rows_in, f"Merge changed the row count from {rows_in:,d} to {len(df):,d}"

    distances = [
        compute_sequence_edit_distances(truth_allele_1, truth_allele_2, tool_allele_1, tool_allele_2)
        for truth_allele_1, truth_allele_2, tool_allele_1, tool_allele_2 in zip(
            df["TruthAlleleSequence: Allele 1"], df["TruthAlleleSequence: Allele 2"],
            df["ToolAlleleSequence: Allele 1"], df["ToolAlleleSequence: Allele 2"])
    ]

    for i, allele in enumerate(("Allele 1", "Allele 2")):
        df[f"SequenceEditDistance: {allele}: {args.tool}"] = [distance[i] for distance in distances]
        # normalize by the truth allele's length, so a 5bp error on a 40bp allele and on a 4kb allele aren't
        # treated as equally wrong. max(1, ...) keeps a zero-length truth allele from dividing by zero.
        df[f"SequenceEditDistanceNormalized: {allele}: {args.tool}"] = [
            None if distance[i] is None else round(distance[i] / max(1, len(sequence_or_none(truth_sequence) or "")), 6)
            for distance, truth_sequence in zip(distances, df[f"TruthAlleleSequence: {allele}"])
        ]

    print_qc(df, args.tool, tool_sequences_df, rows_in)

    # the sequences run to kilobases for VNTRs, so they must not reach the combined variants / alleles tables
    df = df.drop(columns=TRUTH_SEQUENCE_COLUMNS + TOOL_SEQUENCE_COLUMNS
                 + ["ExtractionStatus", "AlleleStatus: Allele 1", "AlleleStatus: Allele 2"])

    output_tsv = args.output_tsv or args.combined_tsv
    df.to_csv(output_tsv, index=False, header=True, sep="\t")
    print(f"Wrote {len(df):,d} rows to {output_tsv}")


if __name__ == "__main__":
    main()
