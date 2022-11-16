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

# paths of other STR catalogs and genomic regions to intersect with
OTHER_STR_CATALOGS = {
    "IlluminaSTRCatalog": "./ref/other/illumina_variant_catalog.sorted.bed.gz",
    "GangSTRCatalog17": "./ref/other/hg38_ver17.fixed.bed.gz",
    "GangSTRCatalog13": "./ref/other/hg38_ver13.fixed.bed.gz",
    "KnownDiseaseAssociatedSTRs": "./ref/other/known_disease_associated_STR_loci.GRCh38.bed.gz",
    "TRFPureRepeats15bp": "./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_15bp.bed.gz",
    "TRFPureRepeats12bp": "./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_12bp.bed.gz",
    "TRFPureRepeats9bp": "./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz",
    "TRFPureRepeats6bp": "./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_6bp.bed.gz",
}

GENOMIC_REGIONS = {
    "SegDupIntervals": "./ref/other/GRCh38GenomicSuperDup.without_decoys.sorted.bed.gz",
}


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
            print(f"Unexpected characters in repeat unit {fields[3]} in line: {line}. Discarding repeat unit...")

    return chrom, start_0based, end_1based, repeat_unit


def create_interval_trees(other_catalogs, counters, show_progress_bar=False, n=None):
    """Create an interval tree for all other catalogs to allow fast overlap checks"""
    all_other_catalog_interval_trees = {}
    for other_catalog_label, is_STR_catalog, other_catalog_path in other_catalogs:
        other_catalog_interval_trees = collections.defaultdict(IntervalTree)
        with gzip.open(os.path.expanduser(other_catalog_path), "rt") as f:

            for i, line in enumerate(f if not show_progress_bar else tqdm.tqdm(f, unit=" bed records")):
                if n is not None and i > n:
                    break
                chrom, start_0based, end_1based, other_catalog_repeat_unit = parse_bed_row(line, is_STR_catalog)
                chrom = chrom.replace("chr", "")
                other_catalog_interval = Interval(start_0based, end_1based, data=other_catalog_repeat_unit)
                other_catalog_interval_trees[chrom].add(other_catalog_interval)
                counters[f"total:{other_catalog_label}"] += 1

        total_count = counters[f"total:{other_catalog_label}"]
        print(f"Parsed {total_count:10,d} intervals from {other_catalog_label:<40s} {other_catalog_path}")

        all_other_catalog_interval_trees[other_catalog_label] = other_catalog_interval_trees

    return all_other_catalog_interval_trees


def process_truth_set_row(truth_set_df, truth_set_row_idx, truth_set_row, other_catalogs, interval_trees, counters, write_beds=False):
    """This method is called for each row in the truth set, and does overlap checks with the interval trees"""
    truth_set_repeat_unit = truth_set_row.Motif
    truth_set_canonical_repeat_unit = compute_canonical_motif(truth_set_repeat_unit, include_reverse_complement=True)
    truth_set_chrom = str(truth_set_row.Chrom).replace("chr", "")
    truth_set_locus_interval = Interval(truth_set_row.Start1Based - 1, truth_set_row.End1Based)
    a1 = int(truth_set_locus_interval.begin)
    a2 = int(truth_set_locus_interval.end)
    counters["total:TruthSetLoci"] += 1
    for other_catalog_label, is_STR_catalog, _ in other_catalogs:
        final_locus_similarity = None
        final_motif_similarity = None
        previous_similarity_rank = 10**6
        if is_STR_catalog:
            locus_similarity_column_name = f"Overlaps{other_catalog_label}: Locus"
        else:
            locus_similarity_column_name = f"Overlaps{other_catalog_label}"
        motif_similarity_column_name = f"Overlaps{other_catalog_label}: Motif"

        # The truth set STR may overlap more than one entry in the other catalog.
        # Look for the entry with the most similarity.
        other_catalog_interval_trees = interval_trees[other_catalog_label]
        overlapping_intervals_from_other_catalog = other_catalog_interval_trees[truth_set_chrom].overlap(truth_set_locus_interval)
        for overlapping_interval_from_other_catalog in overlapping_intervals_from_other_catalog:
            # must overlap by at least 1bp
            if overlapping_interval_from_other_catalog.overlap_size(truth_set_locus_interval) < 1:
                continue

            # compute similarity
            other_catalog_repeat_unit = overlapping_interval_from_other_catalog.data
            if other_catalog_repeat_unit is None:
                # check similarity with interval that is not an STR locus
                motif_similarity = None
                if overlapping_interval_from_other_catalog.contains_interval(truth_set_locus_interval):
                    locus_similarity = "interval contains truth set STR locus"
                    similarity_rank = 0
                elif truth_set_locus_interval.contains_interval(overlapping_interval_from_other_catalog):
                    locus_similarity = "interval contained within truth set STR locus"
                    similarity_rank = 1
                else:
                    locus_similarity = "interval overlaps truth set STR locus"
                    similarity_rank = 2
            else:
                # check similarity with STR from another catalog
                b1 = int(overlapping_interval_from_other_catalog.begin)
                b2 = int(overlapping_interval_from_other_catalog.end)

                other_catalog_canonical_repeat_unit = compute_canonical_motif(other_catalog_repeat_unit, include_reverse_complement=True)
                if a1 == b1 and a2 == b2:
                    locus_similarity = "STR locus has exact same coordinates as truth set STR locus"
                    locus_similarity_rank = 1
                elif overlapping_interval_from_other_catalog.contains_interval(truth_set_locus_interval):
                    locus_similarity = "STR locus contains truth set STR locus"
                    locus_similarity_rank = 2
                elif truth_set_locus_interval.contains_interval(overlapping_interval_from_other_catalog):
                    locus_similarity = "STR locus is contained within truth set STR locus "
                    locus_similarity_rank = 3
                else:
                    locus_similarity = "STR loci overlap"
                    locus_similarity_rank = 4

                if truth_set_repeat_unit == other_catalog_repeat_unit:
                    motif_similarity = "same motif"
                    motif_similarity_rank = 1
                elif truth_set_canonical_repeat_unit == other_catalog_canonical_repeat_unit:
                    motif_similarity = "same motif after normalization"
                    motif_similarity_rank = 10
                else:
                    motif_similarity = "different motif"
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

            if final_locus_similarity is None:
                final_locus_similarity = "no overlap"
                final_motif_similarity = ""

            truth_set_df.at[truth_set_row_idx, locus_similarity_column_name] = final_locus_similarity
            truth_set_df.at[truth_set_row_idx, motif_similarity_column_name] = final_motif_similarity

            counters[f"overlap:{is_STR_catalog}|{other_catalog_label}|{final_locus_similarity}|{final_motif_similarity}"] += 1

            bed_filename = f"STR_{other_catalog_label}_{final_locus_similarity}.bed"
        else:
            truth_set_df.at[truth_set_row_idx, locus_similarity_column_name] = "No" if final_locus_similarity is None else "Yes"
            final_locus_similarity = "no overlap" if final_locus_similarity is None else final_locus_similarity
            final_motif_similarity = ""

            counters[f"overlap:{is_STR_catalog}|{other_catalog_label}|{final_locus_similarity}|"] += 1
            bed_filename = f"region_{other_catalog_label}_{final_locus_similarity}.bed"

        if write_beds and final_locus_similarity != "no overlap" and motif_similarity == "different motif" and any(
                k.lower() in other_catalog_label.lower() for k in ("known", "GangSTR", "Illumina")):
            with open(bed_filename.replace(" ", "_"), "at") as f:
                bed_row_name = f"{truth_set_repeat_unit}:{final_locus_similarity}:{final_motif_similarity}"
                f.write("\t".join(map(str, [truth_set_chrom, a1, a2, bed_row_name, "."])) + "\n")


def print_counters(counters):
    truth_set_total_loci = counters[f"total:TruthSetLoci"]

    def counter_sort_order(x):
        key = x[0]
        count = x[1]
        if key.startswith("total:"):
            return "", count
        key_tokens = key.split("|")

        return key_tokens[1], -len(key_tokens[2]), len(key_tokens[2]), -count

    for key, count in sorted(counters.items(), key=counter_sort_order):
        if key.startswith("total:"):
            continue

        if not key.startswith("overlap:"):
            raise ValueError(f"Unexpected counter key: {key}")

        try:
            is_STR_catalog, other_catalog_label, locus_overlap, motif_overlap = re.sub("^overlap:", "", key).split("|")
        except Exception as e:
            raise ValueError(f"Unexpected counter key: {key}: {e}")

        is_STR_catalog = is_STR_catalog == "True"
        other_catalog_total_loci = counters[f"total:{other_catalog_label}"]

        percent_of_truth_set = f"{100 * count / truth_set_total_loci :5.1f}%" if truth_set_total_loci > 0 else ""
        percent_of_other_catalog = f"{100 * count / other_catalog_total_loci :5.1f}%" if other_catalog_total_loci > 0 else ""

        s = f"{other_catalog_label:30s}: "
        s += f"{count:10,d} out of {truth_set_total_loci:10,d} ({percent_of_truth_set}) truth set loci"
        if locus_overlap == "no overlap":
            s += f": "
        else:
            s += f", and out of {other_catalog_total_loci:10,d} ({percent_of_other_catalog}) {other_catalog_label} loci: "

        s += f"    {locus_overlap:<60s}"
        if is_STR_catalog and motif_overlap:
            print(f"{s}  and has {motif_overlap}")
        else:
            print(s)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", help="Process only the 1st N rows of the truth set TSV. Useful for testing.", type=int)
    p.add_argument("--show-progress-bar", help="Show a progress bar in the terminal when processing variants.", action="store_true")
    p.add_argument("--write-beds", help="Whether to outputu bed files for different overlap categories."
                                       "This is useful for visually checking / troubleshooting results in IGV.", action="store_true")

    p.add_argument("-c", "--catalog", action="append", choices=list(OTHER_STR_CATALOGS) + list(GENOMIC_REGIONS),
                   help="The name of the genomic region or catalog to check overlap with. This option can be specified "
                        "more than once. If not specified, all available catalogs and region lists will be included.")
    p.add_argument("truth_set_tsv_or_bed_path", help="The STR variants TSV file generated by the STR truth set pipeline")
    p.add_argument("output_tsv", nargs="?", help="Optional output file name")
    args = p.parse_args()

    if re.search(".bed(.gz)?", args.truth_set_tsv_or_bed_path):
        truth_set_or_bed_df = pd.read_table(args.truth_set_tsv_or_bed_path, names=[
            "Chrom", "Start0Based", "End1Based", "Motif", "MotifSize"])
        truth_set_or_bed_df.loc[:, "Start1Based"] = truth_set_or_bed_df.Start0Based + 1
    else:
        truth_set_or_bed_df = pd.read_table(args.truth_set_tsv_or_bed_path)


    # filter out non-ref rows
    #len_before = len(truth_set_or_bed_df)
    #truth_set_or_bed_df = truth_set_or_bed_df[truth_set_or_bed_df.End1Based - truth_set_or_bed_df.Start1Based > 0]
    #print(f"Skipped {len_before - len(truth_set_or_bed_df)} loci without matching reference repeats")

    if args.output_tsv:
        output_tsv_path = args.output_tsv
    else:
        output_tsv_path = re.sub(".tsv(.gz)?", "", args.truth_set_tsv_or_bed_path)
        output_tsv_path = re.sub(".bed(.gz)?", "", output_tsv_path)
        output_tsv_path += ".with_overlap_columns.tsv.gz"

    counters = collections.defaultdict(int)

    if args.n:
        truth_set_or_bed_df = truth_set_or_bed_df.iloc[0:args.n, :]

    other_catalogs = [(label, True, path) for label, path in OTHER_STR_CATALOGS.items() if not args.catalog or label in args.catalog]
    other_catalogs += [(label, False, path) for label, path in GENOMIC_REGIONS.items() if not args.catalog or label in args.catalog]

    interval_trees = create_interval_trees(
        other_catalogs, counters, show_progress_bar=args.show_progress_bar, n=args.n)

    # init "Overlaps" columns
    for other_catalog_label, is_STR_catalog, _ in other_catalogs:
        for locus_or_motif in [": Locus", ": Motif"] if is_STR_catalog else [""]:
            truth_set_or_bed_df.at[:, f"Overlaps{other_catalog_label}{locus_or_motif}"] = ""

    # Iterate over all STR variants in the truth set and check overlap
    print(f"Processing rows from {args.truth_set_tsv_or_bed_path}")
    if args.show_progress_bar:
        tqdm.tqdm(truth_set_or_bed_df.iterrows(), total=len(truth_set_or_bed_df), unit=" truth_set records")
    else:
        truth_set_or_bed_df_row_iterator = truth_set_or_bed_df.iterrows()

    for truth_set_row_idx, truth_set_row in truth_set_or_bed_df_row_iterator:
        process_truth_set_row(truth_set_or_bed_df, truth_set_row_idx, truth_set_row,
                              other_catalogs, interval_trees, counters, write_beds=args.write_beds)

    truth_set_or_bed_df.to_csv(output_tsv_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(truth_set_or_bed_df)} rows to the output tsv: {output_tsv_path}")

    print_counters(counters)
    print("Done")


if __name__ == "__main__":
    main()
