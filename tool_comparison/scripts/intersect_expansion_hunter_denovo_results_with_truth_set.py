"""This script takes the raw ExpansionHunterDenovo locus.tsv results table and intersects it with the truth set tsv +
the catalog of all pure reference repeats.
It outputs two TSVs:
truth_set_comparison_table.tsv - one row for every truth set variant, useful for evaluation of EHdn accuracy
EHdn_comparison_table.tsv - one row for every EHdn call, along with a InTruthSet column to indicate which calls matches truth set.
"""

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

import argparse
import collections
import gzip
from intervaltree import Interval, IntervalTree
from pprint import pprint
import os
import pandas as pd
import re
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--syndip-confidence-regions-bed", help="Path of bed file containing the SynDip confidence regions",
                   default="./ref/full.38.bed.gz")
    p.add_argument("--truth-set-variants-tsv", help="Path of truth set variants tsv", default="./STR_truthset.v1.variants.tsv.gz")
    p.add_argument("--reference-repeats-bed", help="Path of bed file containing all STRs in the reference genome",
                   default="./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz")
                    # "./ref/other/repeat_specs_GRCh38_allowing_mismatches.sorted.trimmed.at_least_6bp.bed.gz"
    p.add_argument("-w", "--window-size", type=int, default=500, help="The window size in base pairs. EHdn calls can "
                   "be at most this far away from the true STR locus and still be counted as correct.")
    p.add_argument("expansion_hunter_denovo_results_tsv", nargs="+",
                   help="The output of EHdn ./tool_comparison/results/expansion_hunter_denovo/CHM1_CHM13_WGS2.locus.tsv")
                        #ehdn_tsv = "./tool_comparison/results/expansion_hunter_denovo/CHM1_CHM13_WGS2.locus.tsv"

    args = p.parse_args()
    for path in args.expansion_hunter_denovo_results_tsv + [
        args.syndip_confidence_regions_bed, args.truth_set_variants_tsv, args.reference_repeats_bed
    ]:
        if not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def load_bed_to_interval_trees(bed_path):
    print(f"Loading {bed_path}")
    intervals_trees = collections.defaultdict(IntervalTree)
    fopen = gzip.open if bed_path.endswith("gz") else open
    with fopen(bed_path, "rt") as f:
        for line in tqdm(f, unit=" rows"):
            fields = line.split("\t")
            data = fields[3] if len(fields) > 3 else None
            intervals_trees[fields[0]].add(Interval(int(fields[1]), int(fields[2]), data=data))

    return intervals_trees


def create_function_to_check_EHdn_table_rows_for_overlap(interval_trees, window_size):
    """Returns a function which checks EHdn table rows for whether they overlap any of the intervals defined by the
    given interval trees.

    Args:
         interval_trees (dict): a dictionary that maps chromosome name to IntervalTree.
         window_size (int): the window size in base pairs
    Return:
         function: A function which takes a table row and returns True if that row overlaps any of the intervals.
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
    ehdn_df = pd.read_table(ehdn_results_table_path)
    ehdn_df = ehdn_df[["contig", "start", "end", "motif", "het_str_size"]]  # "num_anc_irrs", "norm_num_anc_irrs",
    ehdn_df = ehdn_df.rename({
        "contig": "Chrom",
        "start": "Start0Based",
        "end": "End1Based",
        "motif": "Motif",
        "het_str_size": "NumRepeats",
    }, axis=1)

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
    ehdn_df.loc[:, "InTruthSet"] = False

    print(f"Read {len(ehdn_df):,d} records from {ehdn_results_table_path}")
    return ehdn_df


def load_truth_set_variants_tsv(truth_set_variants_tsv):
    print(f"Loading {truth_set_variants_tsv}")
    truth_set_df = pd.read_table(truth_set_variants_tsv)
    truth_set_df = truth_set_df[[
        "LocusId", "Locus", "Chrom", "Start1Based", "End1Based",
        "Motif", "CanonicalMotif", "MotifSize",
        "NumRepeatsLongAllele", "RepeatSizeLongAllele (bp)",
        "IsPureRepeat", "IsFoundInReference",

    ]]

    truth_set_df.loc[:, "ReferenceLocusSize (bp)"] = truth_set_df["End1Based"] - truth_set_df["Start1Based"] + 1

    print(f"Read {len(truth_set_df):,d} records from {truth_set_variants_tsv}")
    return truth_set_df


def create_ehdn_output_record(ehdn_call, matches_truth_set_motif=False):
    ehdn_output_record = {}
    for key in "Chrom", "Start1Based", "End1Based", "Motif", "CanonicalMotif", "MotifSize":
        ehdn_output_record[key] = ehdn_call.data[key]

    ehdn_output_record["EHdn NumRepeats"] = ehdn_call.data["NumRepeats"]
    ehdn_output_record["EHdn RepeatSize (bp)"] = ehdn_call.data["RepeatSize (bp)"]
    if ehdn_call.data["MatchingReferenceRepeat"] is None:
        ehdn_output_record["Has Matching Reference Repeat With Same Motif"] = False
    else:
        ehdn_output_record["Has Matching Reference Repeat With Same Motif"] = True
        ehdn_output_record["Matching Reference Repeat Length (bp)"] = ehdn_call.data["MatchingReferenceRepeat"].length()
        ehdn_output_record["Matching Reference Num Repeats"] = ehdn_call.data["MatchingReferenceRepeat"].length() / ehdn_call.data["MotifSize"]

    if matches_truth_set_motif:
        ehdn_output_record["EHdn Concordance With Truth Set"] = "Matches Motif In Truth Set Variant"
    else:
        ehdn_output_record["EHdn Concordance With Truth Set"] = "False Positive"

    return ehdn_output_record


def print_counters(title, counters):
    print(title)
    for key, count in sorted(counters.items(), key=lambda t: -t[1]):
        print(f"  {count:10,d}  {key:10s}")


def main():
    args = parse_args()

    repeats_db_interval_trees = load_bed_to_interval_trees(args.reference_repeats_bed)

    syndip_confidence_interval_trees = load_bed_to_interval_trees(args.syndip_confidence_regions_bed)
    overlaps_syndip_confidence_regions = create_function_to_check_EHdn_table_rows_for_overlap(syndip_confidence_interval_trees, args.window_size)

    truth_set_df = load_truth_set_variants_tsv(args.truth_set_variants_tsv)

    for ehdn_result_table_path in args.expansion_hunter_denovo_results_tsv:
        print("-"*100)

        ehdn_results_df = load_ehdn_results_table(ehdn_result_table_path)
        ehdn_results_df.loc[:, "OverlapsSynDipConfidenceRegions"] = ehdn_results_df.apply(
            overlaps_syndip_confidence_regions, axis=1)
        print(f"Keeping {sum(ehdn_results_df['OverlapsSynDipConfidenceRegions']):,d} out of {len(ehdn_results_df):,d} "
              f"rows that are near syndip confidence intervals")
        ehdn_results_df = ehdn_results_df[ehdn_results_df["OverlapsSynDipConfidenceRegions"]]

        # create interval trees for all EHdn results while also matching them with a nearby reference repeat that
        # has the same canonical motif
        counters = collections.defaultdict(int)
        ehdn_interval_trees = collections.defaultdict(IntervalTree)
        for _, ehdn_results_row in tqdm(ehdn_results_df.iterrows(), total=len(ehdn_results_df)):
            ehdn_metadata = ehdn_results_row.to_dict()

            # find the longest matching reference repeat that has the same canonical motif
            matching_reference_repeat_interval = None
            for reference_repeat_interval in repeats_db_interval_trees[ehdn_results_row.Chrom].overlap(ehdn_results_row.Start1Based - args.window_size, ehdn_results_row.End1Based + args.window_size):
                reference_motif = reference_repeat_interval.data
                if compute_canonical_motif(reference_motif, include_reverse_complement=True) != ehdn_metadata["CanonicalMotif"]:
                    continue

                if matching_reference_repeat_interval is None or matching_reference_repeat_interval.length() < reference_repeat_interval.length():
                    matching_reference_repeat_interval = reference_repeat_interval

            if matching_reference_repeat_interval is not None:
                ehdn_metadata["MatchingReferenceRepeat"] = matching_reference_repeat_interval

                repeat_size_bin = 50 * int(matching_reference_repeat_interval.length() / 50)
                counters[f"Matching reference repeat {repeat_size_bin:03d}bp"] += 1
            else:
                ehdn_metadata["MatchingReferenceRepeat"] = None
                counters[f"No matching reference repeat"] += 1

            ehdn_interval_trees[ehdn_results_row.Chrom].add(
                Interval(ehdn_results_row.Start1Based - args.window_size, ehdn_results_row.End1Based + args.window_size, data=ehdn_metadata)
            )

        print_counters("EHdn matching reference repeats", counters)

        # generate output rows
        ehdn_output_rows = []  # one row per EHdn record
        truth_set_output_rows = []  # one row per truth set record
        counters = collections.defaultdict(int)
        for _, truth_set_row in truth_set_df.iterrows():
            truth_set_output_record = truth_set_row.to_dict()

            truth_set_output_record["EHdn Concordance"] = "No Call"

            # intersect truth set row with EHdn
            overlapping_ehdn_calls = ehdn_interval_trees[truth_set_row.Chrom].overlap(truth_set_row.Start1Based - 1, truth_set_row.End1Based + 1)
            for overlapping_ehdn_call in overlapping_ehdn_calls:
                if len(truth_set_row.CanonicalMotif) == len(overlapping_ehdn_call.data["CanonicalMotif"]):
                    truth_set_output_record["EHdn NumRepeats"] = overlapping_ehdn_call.data["NumRepeats"]
                    truth_set_output_record["EHdn RepeatSize (bp)"] = overlapping_ehdn_call.data["RepeatSize (bp)"]
                    truth_set_output_record["EHdn Concordance"] = "Same Motif Length"

                    if truth_set_row.CanonicalMotif == overlapping_ehdn_call.data["CanonicalMotif"]:
                        repeat_size_bin = 50 * int(truth_set_row["RepeatSizeLongAllele (bp)"] / 50)
                        counters[f"{repeat_size_bin:03d}bp"] += 1

                        overlapping_ehdn_call.data["InTruthSet"] = True
                        truth_set_output_record["EHdn Concordance"] = "Same Motif"

                        ehdn_output_rows.append(
                            create_ehdn_output_record(overlapping_ehdn_call, matches_truth_set_motif=True))
                        break
            truth_set_output_rows.append(truth_set_output_record)
        print_counters("Matching truth set record long allele size (bp)", counters)

        # process the EHdn calls that had no matching truth set variants
        counters = collections.defaultdict(int)
        for ehdn_interval_tree in ehdn_interval_trees.values():
            for ehdn_call in ehdn_interval_tree:
                counters["Total"] += 1
                if ehdn_call.data["InTruthSet"]:
                    counters["TruePositive"] += 1
                    continue

                ehdn_output_rows.append(
                    create_ehdn_output_record(ehdn_call, matches_truth_set_motif=False))
                counters["FalsePositive"] += 1

        print_counters("Truth set concordance", counters)

        # write output tables
        output_path_prefix = re.sub(".gz$", "", ehdn_result_table_path)
        output_path_prefix = re.sub(".tsv$", "", output_path_prefix)
        output_path_prefix = re.sub(".locus$", "", output_path_prefix)

        output_path = f"{output_path_prefix}.truth_set_EHdn_comparison_table.tsv"
        pd.DataFrame(truth_set_output_rows).to_csv(output_path, sep="\t", index=False, header=True)
        print(f"Wrote {len(truth_set_output_rows):,d} rows to {output_path}")

        output_path = f"{output_path_prefix}.EHdn_results_table.tsv"
        pd.DataFrame(ehdn_output_rows).to_csv(output_path, sep="\t", index=False, header=True)
        print(f"Wrote {len(ehdn_output_rows):,d} rows to {output_path}")


if __name__ == "__main__":
    main()