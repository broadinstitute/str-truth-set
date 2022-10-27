"""This script computes statistics about the overlap between the STR truth set vcf and other STR catalogs and
genomic regions of interest.
"""

import argparse
import collections
from intervaltree import Interval, IntervalTree
import gzip
import os
import pandas as pd
import re
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


def parse_bed_row(line, is_STR_catalog):
    fields = line.strip().split("\t")
    chrom = fields[0].strip("chr")
    start_0based = int(fields[1])
    end_1based = int(fields[2])
    repeat_unit = None
    if is_STR_catalog and len(fields) > 3:
        fields[3] = fields[3].replace("(", "").replace(")", "").split("*")[0]
        if re.match("^[ACGTNRYSWKMBDHV]+$", fields[3]):
            repeat_unit = fields[3]
        else:
            print(f"Unexpected characters in repeat unit {fields[3]} in line: {line}. Discarding repeat unit nad setting it to None")

    return chrom, start_0based, end_1based, repeat_unit


def create_interval_trees(other_catalogs, counters, show_progress_bar=False):
    """Create an interval tree for all other catalogs to allow fast overlap checks"""
    all_other_catalog_interval_trees = {}
    for other_catalog_label, is_STR_catalog, other_catalog_path in other_catalogs:
        other_catalog_interval_trees = collections.defaultdict(IntervalTree)
        with gzip.open(os.path.expanduser(other_catalog_path), "rt") as f:

            for line in f if not show_progress_bar else tqdm.tqdm(f, unit=" bed records"):
                chrom, start_0based, end_1based, other_catalog_repeat_unit = parse_bed_row(line, is_STR_catalog)
                chrom = chrom.replace("chr", "")
                other_catalog_interval = Interval(start_0based, end_1based, data=other_catalog_repeat_unit)
                other_catalog_interval_trees[chrom].add(other_catalog_interval)
                counters[f"{other_catalog_label} loci"] += 1

        print(f"Parsed {counters[other_catalog_label+' loci']:,d} intervals from {other_catalog_label}: {other_catalog_path}")
        all_other_catalog_interval_trees[other_catalog_label] = other_catalog_interval_trees

    return all_other_catalog_interval_trees


def process_truth_set_row(truth_set_df, truth_set_row_idx, truth_set_row, other_catalogs, interval_trees, counters):
    """This method is called for each row in the truth set, and performs overlap checks with the given interval trees"""
    truth_set_repeat_unit = truth_set_row.Motif
    truth_set_canonical_repeat_unit = compute_canonical_motif(truth_set_repeat_unit)
    truth_set_chrom = truth_set_row.Chrom.replace("chr", "")
    truth_set_locus_interval = Interval(truth_set_row.Start1Based - 1, truth_set_row.End1Based)
    a1 = int(truth_set_locus_interval.begin)
    a2 = int(truth_set_locus_interval.end)
    counters["Truthset loci"] += 1
    for other_catalog_label, is_STR_catalog, _ in other_catalogs:
        final_locus_similarity = None
        final_motif_similarity = None
        previous_similarity_rank = 10**6
        locus_similarity_column_name = f"Overlaps{other_catalog_label}: Locus" if is_STR_catalog else f"Overlaps{other_catalog_label}"
        motif_similarity_column_name = f"Overlaps{other_catalog_label}: Motif"
        other_catalog_interval_trees = interval_trees[other_catalog_label]
        overlapping_intervals_from_other_catalog = other_catalog_interval_trees[truth_set_chrom].overlap(truth_set_locus_interval)
        for overlapping_interval_from_other_catalog in overlapping_intervals_from_other_catalog:
            # must overlap by more than 1 repeat unit
            if overlapping_interval_from_other_catalog.overlap_size(truth_set_locus_interval) <= len(truth_set_repeat_unit):
                continue

            # compute similarity
            other_catalog_repeat_unit = overlapping_interval_from_other_catalog.data
            if other_catalog_repeat_unit is None:
                motif_similarity = None
                if overlapping_interval_from_other_catalog.contains_interval(truth_set_locus_interval):
                    locus_similarity = "IsWithin"
                    similarity_rank = 0
                elif truth_set_locus_interval.contains_interval(overlapping_interval_from_other_catalog):
                    locus_similarity = "Contains"
                    similarity_rank = 1
                else:
                    locus_similarity = "OverlapsBoundary"
                    similarity_rank = 2
            else:
                b1 = int(overlapping_interval_from_other_catalog.begin)
                b2 = int(overlapping_interval_from_other_catalog.end)

                other_catalog_canonical_repeat_unit = compute_canonical_motif(other_catalog_repeat_unit)
                if a1 == b1 and a2 == b2:
                    locus_similarity = "ExactSameLocus"
                    locus_similarity_rank = 1
                elif overlapping_interval_from_other_catalog.contains_interval(truth_set_locus_interval):
                    locus_similarity = "LocusFromOtherCatalogContainsTruthLocus"
                    locus_similarity_rank = 2
                elif truth_set_locus_interval.contains_interval(overlapping_interval_from_other_catalog):
                    locus_similarity = "TruthLocusContainsLocusFromOtherCatalog"
                    locus_similarity_rank = 3
                else:
                    locus_similarity = "LociOverlap"
                    locus_similarity_rank = 4

                if truth_set_repeat_unit == other_catalog_repeat_unit:
                    motif_similarity = "ExactSameMotif"
                    motif_similarity_rank = 1
                elif truth_set_canonical_repeat_unit == other_catalog_canonical_repeat_unit:
                    motif_similarity = "ShiftedMotif"
                    motif_similarity_rank = 10
                else:
                    motif_similarity = "CompletelyDifferentMotif"
                    motif_similarity_rank = 1000

                similarity_rank = motif_similarity_rank * locus_similarity_rank - 1

            if similarity_rank < previous_similarity_rank:
                previous_similarity_rank = similarity_rank
                final_motif_similarity = motif_similarity
                final_locus_similarity = locus_similarity
                if similarity_rank == 0:
                    # can't do any better than this
                    break

        if is_STR_catalog:
            if final_motif_similarity is None:
                final_motif_similarity = "NoOverlappingLoci"

            if final_locus_similarity is None:
                final_locus_similarity = "NoOverlappingLoci"
        else:
            if final_locus_similarity is None:
                final_locus_similarity = "NoOverlap"

        truth_set_df.at[truth_set_row_idx, locus_similarity_column_name] = final_locus_similarity
        if is_STR_catalog:
            truth_set_df.at[truth_set_row_idx, motif_similarity_column_name] = final_motif_similarity
            counters[f"{other_catalog_label} vs Truthset: {final_motif_similarity}|{final_locus_similarity}"] += 1
        else:
            counters[f"{other_catalog_label} vs Truthset: {final_locus_similarity}"] += 1


def print_counters(counters):
    def counter_sort_order(x):
        split_by_sep1 = x[0].split(":")
        #split_by_sep2 = split_by_sep1[1].split("|") if len(split_by_sep1) > 1 else "X"
        return (-1 * x[0].endswith(" loci"), split_by_sep1[0], -x[1]) #, split_by_sep2[1] if len(split_by_sep2) > 1 else "X")

    total_truth_set_loci = counters[f"Truthset loci"]
    for key, value in sorted(counters.items(), key=counter_sort_order):
        percent_of_truth_set = f"{100 * value / total_truth_set_loci :5.1f}%"

        if key.endswith("loci"):
            print(f"{value:10,d} {key}")
        else:
            s = f"{value:10,d}  {key:100s}   ({percent_of_truth_set} of {total_truth_set_loci:10,d} truth_set loci)"
            if " vs Truthset" in key and not key.endswith(": NoOverlap|None"):
                other_catalog_label = key.split(" vs Truthset")[0]
                total_loci_in_other_catalog = counters[f"{other_catalog_label} loci"]
                percent_of_other_catalog = f"{100 * value / total_loci_in_other_catalog :5.1f}%"
                s += f"    ({percent_of_other_catalog} of {total_loci_in_other_catalog:10,d} {other_catalog_label} loci)"
            print(s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", help="Process only the 1st N rows of the truth set TSV. Useful for testing.", type=int)
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    p.add_argument("truth_set_tsv", help="The STR variants TSV file generated by the STR truth set pipeline" )
    p.add_argument("output_tsv", nargs="?", help="Optional output file name")
    args = p.parse_args()

    truth_set_df = pd.read_table(args.truth_set_tsv)

    # filter out non-ref rows
    len_before = len(truth_set_df)
    truth_set_df = truth_set_df[truth_set_df.End1Based - truth_set_df.Start1Based > 0]
    print(f"Skipped {len_before - len(truth_set_df)} loci without matching reference repeats")

    output_tsv_path = args.output_tsv if args.output_tsv else args.truth_set_tsv.replace(".tsv", "") + ".with_overlap_columns.tsv"

    counters = collections.defaultdict(int)

    # load other catalogs into IntervalTrees
    other_catalogs = [
        ("KnownPathogenicSTRCatalog", True, "./ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz"),
        ("IlluminaSTRCatalog", True, "./ref/other/illumina_variant_catalog.sorted.bed.gz"),
        ("GangSTRCatalog", True, "./ref/other/hg38_ver17.fixed.bed.gz"),
        ("SegDupIntervals", False, "./ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz"),
        ("ExonsFromGencode40", False, "./ref/other/gencode.v40.exons.bed.gz"),
        ("CDSFromGencode40", False, "./ref/other/gencode.v40.CDS.bed.gz"),
        ("ExonsFromMANE", False, "./ref/other/MANE.v1.exons.bed.gz"),
        ("CDSFromMANE", False, "./ref/other/MANE.v1.CDS.bed.gz"),
    ]

    if args.n:
        truth_set_df = truth_set_df.iloc[0:args.n, :]
        #other_catalogs = other_catalogs[:1]

    interval_trees = create_interval_trees(
        other_catalogs, counters, show_progress_bar=args.show_progress_bar)

    # init "Overlaps" columns
    for other_catalog_label, is_STR_catalog, _ in other_catalogs:
        for locus_or_motif in [": Locus", ": Motif"] if is_STR_catalog else [""]:
            truth_set_df.at[:, f"Overlaps{other_catalog_label}{locus_or_motif}"] = ""

    # Iterate over all STR variants in the truth set and check overlap
    print(f"Processing rows from {args.truth_set_tsv}")
    if args.show_progress_bar:
        tqdm.tqdm(truth_set_df.iterrows(), total=len(truth_set_df), unit=" truth_set records")
    else:
        truth_set_df_row_iterator = truth_set_df.iterrows()

    for truth_set_row_idx, truth_set_row in truth_set_df_row_iterator:
        process_truth_set_row(
            truth_set_df, truth_set_row_idx, truth_set_row, other_catalogs, interval_trees, counters)

    truth_set_df.to_csv(output_tsv_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(truth_set_df)} rows to the output tsv: {output_tsv_path}")

    print_counters(counters)
    print("Done")


if __name__ == "__main__":
    main()

#%%