"""This script takes the raw ExpansionHunterDenovo locus.tsv results table and intersects it with the truth set tsv
as well as the catalog of all pure STRs in the reference genome.

It then outputs two TSVs:

{prefix}.truth_set_EHdn_comparison_table.tsv - one row for every truth set allele, for evaluation of EHdn accuracy.
{prefix}.EHdn_results_table.tsv - one row for every EHdn call, along with a new "InTruthSet" column to indicate 
    which calls match a truth set allele. This is useful for evaluating false-positive calls.
"""

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

import argparse
import collections
import gzip
from intervaltree import Interval, IntervalTree
import os
import pandas as pd
import re
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--syndip-confidence-regions-bed", default="./ref/full.38.bed.gz",
                   help="Path of bed file containing the SynDip confidence regions, used to filter out EHdn results in "
                        "regions that are excluded from the STR truth set")
    p.add_argument("--truth-set-alleles-tsv", default="./STR_truth_set.v1.alleles.tsv.gz",
                   help="Path of truth set alleles tsv")
    p.add_argument("--reference-repeats-bed", help="Path of bed file containing all STRs in the reference genome",
                   default="./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz")
                    # "./ref/other/repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_6bp.bed.gz"
    p.add_argument("-w", "--window-size", type=int, default=600, help="The window size in base pairs. EHdn calls can "
                   "be at most this far away from the true STR locus and still be counted as correct.")
    p.add_argument("expansion_hunter_denovo_results_tsv", nargs="+",
                   help="The output of EHdn ./tool_comparison/results/expansion_hunter_denovo/CHM1_CHM13_WGS2.locus.tsv")
                        #ehdn_tsv = "./tool_comparison/results/expansion_hunter_denovo/CHM1_CHM13_WGS2.locus.tsv"
    args = p.parse_args()

    # validate file path args
    for path in args.expansion_hunter_denovo_results_tsv + [
        args.syndip_confidence_regions_bed, args.truth_set_variants_tsv, args.reference_repeats_bed
    ]:
        if not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def load_bed_to_interval_trees(bed_path):
    """Load a bed file and return a dictionary of IntervalTrees - one for each chromosome."""
    print(f"Loading {bed_path}")
    intervals_trees = collections.defaultdict(IntervalTree)
    fopen = gzip.open if bed_path.endswith("gz") else open
    with fopen(bed_path, "rt") as f:
        for line in tqdm(f, unit=" rows", unit_scale=True):
            fields = line.split("\t")
            chrom = fields[0]
            start_0based = int(fields[1])
            end_1based = int(fields[2])
            data = fields[3] if len(fields) > 3 else None
            intervals_trees[chrom].add(Interval(start_0based, end_1based, data=data))

    return intervals_trees


def create_function_to_check_ehdn_table_rows_for_overlap_with_interval_trees(interval_trees, window_size):
    """Returns a function which checks EHdn table rows for whether they overlap any of the intervals defined by the
    given interval trees and after widening the EHdn interval by 'window_size' base pairs in each direction.

    Args:
         interval_trees (dict): a dictionary that maps chromosome name to IntervalTree.
         window_size (int): the window size in base pairs by which to widen the EHdn interval in each direction.
            This way, the EHdn interval can be at most this far away from an interval in the IntervalTrees and still be
            considered as overlapping.
    Return:
         function: A function which takes a table row and returns True if that row is closer than window_size base pairs
            to one of the intervals.
    """
    def overlaps_interval_trees(ehdn_results_table_row):
        """Returns true if the given EHdn results table row overlaps any of the intervals."""

        if ehdn_results_table_row.Chrom not in interval_trees:
            raise ValueError(f"Unexpected chromosome: {ehdn_results_table_row.Chrom}")

        if interval_trees[ehdn_results_table_row.Chrom].overlaps(
                ehdn_results_table_row.Start0Based - window_size,
                ehdn_results_table_row.End1Based + window_size):
            return True

        return False

    return overlaps_interval_trees


def load_ehdn_results_table(ehdn_results_table_path):
    """Loads the combined EHdn results locus table. Then filters rows to leave only the autosomes and chrX.
    Then, computes additional columns and returns the table.
    """
    ehdn_df = pd.read_table(ehdn_results_table_path)
    ehdn_df = ehdn_df[["contig", "start", "end", "motif", "num_anc_irrs", "norm_num_anc_irrs", "het_str_size"]]
    ehdn_df = ehdn_df.rename({
        "contig": "Chrom",
        "start": "Start0Based",
        "end": "End1Based",
        "motif": "Motif",
        "num_anc_irrs": "NumAnchoredIRRs",
        "norm_num_anc_irrs": "NumAnchoredIRRsNormalized",
        "het_str_size": "NumRepeats",
    }, axis=1)

    assert all(ehdn_df.Start0Based.astype(int) < ehdn_df.End1Based.astype(int))

    ehdn_df = ehdn_df[
        ~ehdn_df.Chrom.str.startswith("chrUn_") &
        ~ehdn_df.Chrom.str.endswith("_random") &
        ~ehdn_df.Chrom.str.endswith("_alt") &
        ~ehdn_df.Chrom.str.endswith("chrM") &
        ~ehdn_df.Chrom.str.endswith("chrY")
    ]
    ehdn_df.loc[:, "Start1Based"] = ehdn_df.Start0Based + 1
    ehdn_df.loc[:, "CanonicalMotif"] = ehdn_df.Motif.apply(lambda s: compute_canonical_motif(s, include_reverse_complement=True))
    ehdn_df.loc[:, "MotifSize"] = ehdn_df.Motif.str.len()
    ehdn_df.loc[:, "RepeatSize (bp)"] = ehdn_df.NumRepeats * ehdn_df.MotifSize
    ehdn_df.loc[:, "ReferenceRegion"] = ehdn_df.Chrom + ":" + ehdn_df.Start1Based.astype(str) + "-" + ehdn_df.End1Based.astype(str)

    print(f"Read {len(ehdn_df):,d} records from {ehdn_results_table_path}")
    return ehdn_df


def load_truth_set_variants_tsv(truth_set_variants_tsv):
    """Loads and returns the truth set tsv after discarding irrelevant columns and adding a "ReferenceLocusSize (bp)"
    column.
    """

    print(f"Loading {truth_set_variants_tsv}")
    truth_set_variants_df = pd.read_table(truth_set_variants_tsv)
    truth_set_variants_df = truth_set_variants_df[[
        "LocusId", "Locus", "Chrom", "Start1Based", "End1Based",
        "Motif", "CanonicalMotif", "MotifSize", "INS_or_DEL",
        "NumRepeatsInReference", "NumRepeats", "RepeatSize (bp)",
        "IsPureRepeat", "IsFoundInReference", "SummaryString",
    ]]

    # filter to expansions and motif sizes that are smaller than the MAX_REPEAT_UNIT_LENGTH set in the EHdn pipeline
    truth_set_variants_df = truth_set_variants_df[truth_set_variants_df["INS_or_DEL"].str.contains("INS")]
    truth_set_variants_df = truth_set_variants_df[truth_set_variants_df["MotifSize"] <= 50]

    truth_set_variants_df.loc[:, "ReferenceLocusSize (bp)"] = truth_set_variants_df["End1Based"] - truth_set_variants_df["Start1Based"] + 1

    print(f"Read {len(truth_set_variants_df):,d} STR expansion records from {truth_set_variants_tsv}")
    return truth_set_variants_df


def create_truth_set_output_record(truth_set_row, matching_ehdn_calls):
    truth_set_output_record = truth_set_row.to_dict()

    if len(matching_ehdn_calls) == 0:
        truth_set_output_record["EHdn Concordance"] = "No Call"
        truth_set_output_record["EHdn ReferenceRegion"] = None

    elif len(matching_ehdn_calls) == 1:
        truth_set_output_record["EHdn Concordance"] = f"Matching EHdn Call"
        truth_set_output_record["EHdn ReferenceRegion"] = matching_ehdn_calls[0].data["ReferenceRegion"]

    elif len(matching_ehdn_calls) > 1:
        # just to see if there are any truth set loci that are near multiple EHdn calls with the same motif
        truth_set_output_record["EHdn Concordance"] = f"Matching EHdn Call"

        # find the closest EHdn region
        print(f"Truth set locus", truth_set_output_record["SummaryString"], f"has {len(matching_ehdn_calls)} matching EHdn calls:")
        min_distance = 10**9
        for matching_ehdn_call in matching_ehdn_calls:
            a = int(truth_set_output_record["Start1Based"]) - int(matching_ehdn_call.data["End1Based"])
            b = int(matching_ehdn_call.data["Start0Based"]) - int(truth_set_output_record["End1Based"])

            if a > 0:    # truth set locus is to the right of the EHdn locus
                print(" "*10, f"{a:,d}bp distance between EHdn call @", matching_ehdn_call.data["ReferenceRegion"],
                      "and truth set locus", truth_set_output_record["Locus"])
                if a < min_distance:
                    min_distance = a
                    truth_set_output_record["EHdn ReferenceRegion"] = matching_ehdn_call.data["ReferenceRegion"]

            elif b > 0:  # EHdn locus is to the right of the truth set locus
                print(" "*10, f"{b:,d}bp distance between truth set locus", truth_set_output_record["Locus"],
                      "and EHdn call @", matching_ehdn_call.data["ReferenceRegion"])
                if b < min_distance:
                    min_distance = b
                    truth_set_output_record["EHdn ReferenceRegion"] = matching_ehdn_call.data["ReferenceRegion"]

            else:  # loci overlap
                print(" "*10, "0bp distance between truth set locus", truth_set_output_record["Locus"],
                      "and EHdn call @", matching_ehdn_call.data["ReferenceRegion"], "since they overlap")
                min_distance = 0
                truth_set_output_record["EHdn ReferenceRegion"] = matching_ehdn_call.data["ReferenceRegion"]

    return truth_set_output_record


def create_ehdn_output_record(ehdn_call):
    """Computes output table rows"""

    ehdn_output_record = {}
    for key in "Chrom", "Start1Based", "End1Based", "Motif", "CanonicalMotif", "MotifSize":
        ehdn_output_record[key] = ehdn_call.data[key]

    ehdn_output_record["EHdn Locus"] = "".join(map(str, [
        ehdn_call.data["Chrom"], ":", ehdn_call.data["Start1Based"], "-", ehdn_call.data["End1Based"]
    ]))
    ehdn_output_record["EHdn LocusSize (bp)"] = int(ehdn_call.data["End1Based"]) - int(ehdn_call.data["Start0Based"])
    ehdn_output_record["EHdn NumRepeats"] = ehdn_call.data["NumRepeats"]
    ehdn_output_record["EHdn RepeatSize (bp)"] = ehdn_call.data["RepeatSize (bp)"]
    ehdn_output_record["EHdn NumAnchoredIRRs"] = ehdn_call.data["NumAnchoredIRRs"]
    ehdn_output_record["EHdn NumAnchoredIRRsNormalized"] = ehdn_call.data["NumAnchoredIRRsNormalized"]

    ehdn_output_record["HasMatchingReferenceLocus"] = ehdn_call.data["MatchingReferenceLocus"] is not None
    ehdn_output_record["MatchingReferenceLocus"] = ehdn_call.data["MatchingReferenceLocus"]
    ehdn_output_record["MatchingReferenceLocusSize (bp)"] = ehdn_call.data["MatchingReferenceLocusSize (bp)"]

    matching_truth_set_row = ehdn_call.data["MatchingTruthSetRow"]
    if matching_truth_set_row is not None:
        ehdn_output_record["EHdn Concordance With Truth Set"] = "Matching Expansion In Truth Set"
        ehdn_output_record["TruthSet NumRepeats"] = matching_truth_set_row["NumRepeatsLongAllele"]
        ehdn_output_record["TruthSet RepeatSize (bp)"] = matching_truth_set_row["RepeatSizeLongAllele (bp)"]
        ehdn_output_record["TruthSet IsPureRepeat"] = matching_truth_set_row["IsPureRepeat"]
    else:
        ehdn_output_record["EHdn Concordance With Truth Set"] = "False Positive"
        ehdn_output_record["TruthSet NumRepeats"] = None
        ehdn_output_record["TruthSet RepeatSize (bp)"] = None
        ehdn_output_record["TruthSet IsPureRepeat"] = None

    return ehdn_output_record


def print_counters(title, counters):
    """Utility method to print a dictionary of countes"""
    print(title)
    total = counters.get("Total")
    for key, count in sorted(counters.items(), key=lambda t: -t[1]):
        if total is None:
            print(f"  {count:10,d}  {key:10s}")
        elif key != "Total":
            print(f"  {count:10,d}  ({100*count/total:0.1f}%)  {key:10s}")


def find_longest_matching_reference_repeat(ehdn_call, reference_repeats_interval_trees, window_size):
    """Utility method that looks for a reference genome repeat locus that is near an EHdn call and has the same motif"""

    # find the longest overlapping reference repeat that has the same canonical motif
    matching_reference_repeat_interval = None
    for reference_repeat_interval in reference_repeats_interval_trees[ehdn_call.data["Chrom"]].overlap(
            ehdn_call.data["Start1Based"] - window_size, ehdn_call.data["End1Based"] + window_size):

        reference_motif = str(reference_repeat_interval.data)
        if (matching_reference_repeat_interval is None
            or matching_reference_repeat_interval.length() < reference_repeat_interval.length()) and (
                compute_canonical_motif(reference_motif, include_reverse_complement=True) == ehdn_call.data["CanonicalMotif"]
        ):
            matching_reference_repeat_interval = reference_repeat_interval

    if matching_reference_repeat_interval is None:
        return None, 0
    else:
        chrom = ehdn_call.data["Chrom"]
        start_0based = matching_reference_repeat_interval.begin
        end_1based = matching_reference_repeat_interval.end

        return (
            f"{chrom}:{start_0based}-{end_1based}",
            matching_reference_repeat_interval.end - matching_reference_repeat_interval.begin
        )


def generate_output_table(
        ehdn_result_table_path,
        truth_set_df,
        syndip_confidence_interval_trees,
        reference_repeats_interval_trees,
        window_size):

    print("-"*100)
    overlaps_syndip_confidence_regions = create_function_to_check_ehdn_table_rows_for_overlap_with_interval_trees(
        syndip_confidence_interval_trees, window_size)

    # load ehdn_results_df
    print(f"Loading {ehdn_result_table_path}")
    ehdn_results_df = load_ehdn_results_table(ehdn_result_table_path)

    # filter ehdn_results_df to rows that are close to syndip confidence intervals
    print(f"Filtering {ehdn_result_table_path} to syndip confidence regions")
    ehdn_results_df.loc[:, "OverlapsSynDipConfidenceRegions"] = ehdn_results_df.apply(
        overlaps_syndip_confidence_regions, axis=1)
    num_overlap_syndip_confidence_regions = sum(ehdn_results_df['OverlapsSynDipConfidenceRegions'])
    print(f"Keeping {num_overlap_syndip_confidence_regions:,d} out of {len(ehdn_results_df):,d} "
          f"({100*num_overlap_syndip_confidence_regions/len(ehdn_results_df):0.1f}%) rows that are near "
          f"syndip confidence intervals")
    ehdn_results_df = ehdn_results_df[ehdn_results_df["OverlapsSynDipConfidenceRegions"]]

    # create interval trees for all EHdn calls
    print("Create interval trees for EHdn calls")
    ehdn_interval_trees = collections.defaultdict(IntervalTree)
    for _, ehdn_results_row in tqdm(ehdn_results_df.iterrows(), total=len(ehdn_results_df), unit=" EHdn calls"):
        ehdn_results_row_dict = ehdn_results_row.to_dict()
        ehdn_results_row_dict["MatchingTruthSetRows"] = []  # used to keep track of any adjacent truth set rows with matching motifs
        ehdn_interval = Interval(
            ehdn_results_row.Start1Based - 1, ehdn_results_row.End1Based + 1, data=ehdn_results_row_dict)
        ehdn_interval_trees[ehdn_results_row.Chrom].add(ehdn_interval)

    # iterate over all truth set rows and match them to EHdn calls while generating truth_set_ehdn_comparison_table_rows
    print("Match truth set alleles with EHdn calls")
    truth_set_ehdn_comparison_table_rows = []  # one row per truth set record
    for _, truth_set_row in tqdm(truth_set_df.iterrows(), total=len(truth_set_df), unit=" truth set alleles"):

        # intersect truth set row with EHdn
        matching_ehdn_calls = []
        overlapping_ehdn_calls = ehdn_interval_trees[truth_set_row.Chrom].overlap(
            truth_set_row.Start1Based - window_size, truth_set_row.End1Based + window_size)
        for overlapping_ehdn_call in overlapping_ehdn_calls:
            if truth_set_row.CanonicalMotif == overlapping_ehdn_call.data["CanonicalMotif"]:
                matching_ehdn_calls.append(overlapping_ehdn_call)
                overlapping_ehdn_call.data["MatchingTruthSetRows"].append(truth_set_row)

        r = create_truth_set_output_record(truth_set_row, matching_ehdn_calls)
        truth_set_ehdn_comparison_table_rows.append(r)

    # generate ehdn_results_table_rows
    print("Generate EHdn results table")
    counters = collections.defaultdict(int)
    ehdn_results_table_rows = []               # one row per EHdn record
    ehdn_call_generator = (
        ehdn_call for ehdn_interval_tree in ehdn_interval_trees.values() for ehdn_call in ehdn_interval_tree)

    for ehdn_call in tqdm(ehdn_call_generator, unit=" EHdn call", total=len(ehdn_results_df)):
        counters["Total"] += 1

        if len(ehdn_call.data["MatchingTruthSetRows"]) > 0:

            # get the reference locus from the overlapping truth set allele with the largest expansion
            matching_truth_set_row_with_largest_expansion = None
            for matching_truth_set_row in ehdn_call.data["MatchingTruthSetRows"]:
                if matching_truth_set_row_with_largest_expansion is None or (
                        matching_truth_set_row["NumRepeatsLongAllele"] >
                        matching_truth_set_row_with_largest_expansion["NumRepeatsLongAllele"]):
                    matching_truth_set_row_with_largest_expansion = matching_truth_set_row

            ehdn_call.data["MatchingTruthSetRow"] = matching_truth_set_row_with_largest_expansion
            ehdn_call.data["MatchingReferenceLocus"] = matching_truth_set_row_with_largest_expansion["Locus"]
            ehdn_call.data["MatchingReferenceLocusSize (bp)"] = (
                    matching_truth_set_row_with_largest_expansion["NumRepeatsInReference"] * len(
                    matching_truth_set_row_with_largest_expansion["Motif"])
            )

            counters["TruePositive"] += 1

        else:
            ehdn_call.data["MatchingTruthSetRow"] = None
            ehdn_call.data["MatchingReferenceLocus"], ehdn_call.data["MatchingReferenceLocusSize (bp)"] = find_longest_matching_reference_repeat(
                ehdn_call, reference_repeats_interval_trees, window_size)

            if ehdn_call.data["MatchingReferenceLocus"]:
                counters["False Positive With Matching Reference Locus"] += 1
            else:
                counters["False Positive Without Matching Reference Locus"] += 1

        output_row = create_ehdn_output_record(ehdn_call)
        ehdn_results_table_rows.append(output_row)

    # Matching them with a nearby reference repeat that has the same canonical motif
    print_counters("EHdn vs truth set repeats", counters)

    # write output tables
    print("Write output tables")
    output_path_prefix = re.sub(".gz$", "", ehdn_result_table_path)
    output_path_prefix = re.sub(".tsv$", "", output_path_prefix)
    output_path_prefix = re.sub(".locus$", "", output_path_prefix)

    output_path = f"{output_path_prefix}.truth_set_EHdn_comparison_table.tsv"
    pd.DataFrame(truth_set_ehdn_comparison_table_rows).to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(truth_set_ehdn_comparison_table_rows):,d} rows to {output_path}")

    output_path = f"{output_path_prefix}.EHdn_results_table.with_truth_set_concordance.tsv"
    pd.DataFrame(ehdn_results_table_rows).to_csv(output_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(ehdn_results_table_rows):,d} rows to {output_path}")


def main():
    args = parse_args()

    reference_repeats_interval_trees = load_bed_to_interval_trees(args.reference_repeats_bed)

    syndip_confidence_interval_trees = load_bed_to_interval_trees(args.syndip_confidence_regions_bed)

    truth_set_df = load_truth_set_variants_tsv(args.truth_set_variants_tsv)

    for ehdn_result_table_path in sorted(args.expansion_hunter_denovo_results_tsv):
        generate_output_table(
            ehdn_result_table_path,
            truth_set_df,
            syndip_confidence_interval_trees,
            reference_repeats_interval_trees,
            args.window_size)


if __name__ == "__main__":
    main()