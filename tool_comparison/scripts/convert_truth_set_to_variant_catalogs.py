"""
Create ExpansionHunter and GangSTR variant catalogs for all loci in the STR truth set as well as an equal
number of loci from other places in the genome to use as true-negatives.
"""

import argparse
import collections
from intervaltree import Interval, IntervalTree
import json
import os
import pandas as pd
from pprint import pprint
import pybedtools
import random
import re
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

random.seed(1)

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--expansion-hunter-loci-per-run", type=int, default=500, help="ExpansionHunter batch size. "
                   "The set of all STR loci in the truth set will be split into variant catalogs of this size.")
    p.add_argument("--gangstr-loci-per-run", type=int, default=10000, help="GangSTR batch size. "
                   "The set of all STR loci in the truth set will be split into repeat spec files of this size.")
    p.add_argument("--high-confidence-regions-bed", default="./ref/full.38.bed.gz",
                    help="Path of the SynDip high-confidence regions .bed file")
    p.add_argument("--all-repeats-bed", default="./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz",
                   help="Path of bed file containing all repeats in the reference genome generated using a tool like "
                        "TandemRepeatFinder")
    p.add_argument("--skip-negative-loci", action="store_true", help="Skip generating negative loci")
    p.add_argument("--skip-eh-catalog", action="store_true", help="Skip generating an ExpansionHunter catalog")
    p.add_argument("--skip-gangstr-catalog", action="store_true", help="Skip generating a GangSTR catalog")
    p.add_argument("--skip-hipstr-catalog", action="store_true", help="Skip generating a HipSTR catalog")
    p.add_argument("--output-dir", default="./tool_comparison/variant_catalogs/", help="Directory where to write output files")
    p.add_argument("truth_set_variants_tsv_or_bed_path", nargs="?", default="STR_truth_set.v1.variants.tsv",
                   help="Path of the STR truth set .variants.tsv or of an arbitrary bed file")
    args = p.parse_args()

    for path in args.truth_set_variants_tsv_or_bed_path, args.all_repeats_bed, args.high_confidence_regions_bed:
        if not os.path.isfile(path):
            p.error(f"{path} not found")

    return args


def compute_row_key(motif):
    # above this threshold, stratify loci only by motif size and not hte motif itself
    if len(motif) == 2:
        key = (len(motif), motif)
    elif len(motif) == 3 or len(motif) == 4:
        canonical_motif = compute_canonical_motif(motif, include_reverse_complement=True)
        key = (len(canonical_motif), canonical_motif)
    elif len(motif) > 30:
        key = (30, None)
    else:
        key = (len(motif), None)

    return key


def generate_set_of_positive_loci(truth_set_df):
    """Convert the truth set table into a set of positive loci, as well as a dictionary of counters stratified by motif,
    and interval trees for efficiently computing interval overlaps with loci in the truth set.
    """

    positive_loci = set()  # loci with STR variants in CHM1-CHM13
    positive_loci_counters = collections.defaultdict(int)
    truth_set_loci_interval_trees = collections.defaultdict(IntervalTree)
    for _, row in truth_set_df.iterrows():
        motif = row.Motif
        chrom = row.Chrom.replace("chr", "")
        start_0based = int(row.Start1Based) - 1
        end = int(row.End1Based)
        end -= (end - start_0based) % len(motif)
        if end != int(row.End1Based):
            print("Trimmed", row)

        key = compute_row_key(motif)
        positive_loci.add((chrom, start_0based, end, motif))
        positive_loci_counters[key] += 1

        truth_set_loci_interval_trees[chrom].add(Interval(start_0based, end, data=motif))

    return positive_loci, positive_loci_counters, truth_set_loci_interval_trees


def load_all_pure_repeats_in_syndip_confidence_regions(all_repeats_bed_path, high_confidence_repeats_bed_path):
    """Load the bed file containing all pure STR repeats in the reference genome (defined by running a tool like
    TandemRepeatFinder) and filter it to SynDip high-confidence regions.
    This will be used to define the set of non-variable (ie. negative regions)
    """

    # get all pure repeats in GRCh38 that are within syndip high confidence regions
    pure_repeats_bedtool = pybedtools.BedTool(os.path.expanduser(all_repeats_bed_path))
    all_pure_repeats_in_syndip_high_confidence_regions_iter = pure_repeats_bedtool.intersect(
        os.path.expanduser(high_confidence_repeats_bed_path), f=1, wa=True, u=True)

    all_pure_repeats_in_syndip_high_confidence_regions = []
    for row in all_pure_repeats_in_syndip_high_confidence_regions_iter:
        motif = row.name
        chrom = row.chrom.replace("chr", "")
        start_0based = int(row.start)
        end = int(row.end)
        end -= (end - start_0based) % len(motif)
        if int(end) != int(row.end):
            print("Trimmed", row)

        all_pure_repeats_in_syndip_high_confidence_regions.append((chrom, start_0based, end, motif))

    return all_pure_repeats_in_syndip_high_confidence_regions


def generate_set_of_negative_loci(
        all_pure_repeats_in_syndip_high_confidence_regions, truth_set_loci_interval_trees, positive_loci_counters):

    # shuffle so that negative examples are distributed throughout the genome
    random.shuffle(all_pure_repeats_in_syndip_high_confidence_regions)

    # compute negative loci
    negative_loci = set()
    negative_loci_counters = collections.defaultdict(int)
    enough_negative_loci = collections.defaultdict(bool)
    for chrom, start_0based, end, motif in all_pure_repeats_in_syndip_high_confidence_regions:
        key = compute_row_key(motif)
        if enough_negative_loci[key]:
            continue

        # skip variant loci (ie. loci that overlap the truth set)
        if truth_set_loci_interval_trees[chrom].overlap(Interval(start_0based, end, data=motif)):
            continue

        negative_loci.add((chrom, start_0based, end, motif))
        negative_loci_counters[key] += 1
        if negative_loci_counters[key] >= positive_loci_counters[key]:
            enough_negative_loci[key] = True

            # check if finished
            if all(enough_negative_loci.values()):
                break

    return negative_loci


def write_expansion_hunter_variant_catalogs(locus_set, output_path_prefix, loci_per_run=None):
    variant_catalog = []
    for chrom, start_0based, end_1based, motif in locus_set:
        variant_catalog.append({
            "LocusId": f"{chrom}-{start_0based}-{end_1based}-{motif}",
            "ReferenceRegion": f"{chrom}:{start_0based}-{end_1based}",
            "LocusStructure": f"({motif})*",
            "VariantType": "Repeat",
        })

    if loci_per_run is None:
        batches = [variant_catalog]
    else:
        batches = [
            variant_catalog[i:i+loci_per_run] for i in range(0, len(variant_catalog), loci_per_run)
        ]
    for batch_i, current_variant_catalog in enumerate(batches):
        with open(f"{output_path_prefix}.{batch_i+1:03d}_of_{len(batches):03d}.json", "wt") as f:
            json.dump(current_variant_catalog, f, indent=3)

    print(f"Wrote {len(batches):,d} ExpansionHunter variant catalogs to {output_path_prefix}*.json")


def write_gangstr_or_hipstr_repeat_specs(locus_set, output_path_prefix, gangstr=False, hipstr=False, loci_per_run=None):
    locus_list = list(locus_set)

    if (not gangstr and not hipstr) or (gangstr and hipstr):
        raise ValueError(f"Can not set gangstr={gangstr} and hipstr={hipstr}")

    # write out GangSTR repeat specs
    if loci_per_run is None:
        batches = [locus_list]
    else:
        batches = [
            locus_list[i:i+loci_per_run] for i in range(0, len(locus_list), loci_per_run)
        ]

    for batch_i, current_repeat_specs in enumerate(batches):
        with open(f"{output_path_prefix}.{batch_i+1:03d}_of_{len(batches):03d}.bed", "wt") as f:
            for chrom, start_0based, end_1based, motif in current_repeat_specs:
                chrom = f"chr" + chrom.replace("chr", "")
                if gangstr:
                    output_fields = [chrom, start_0based + 1, end_1based, len(motif), motif]
                elif hipstr:
                    if len(motif) > 9 or (end_1based - start_0based) <= 1:
                        # HipSTR doesn't support motifs longer than 9bp or loci where start_1based == stop_1based
                        # (https://github.com/tfwillems/HipSTR/blob/master/src/region.cpp#L33-L35)
                        continue

                    output_fields = [
                        chrom, start_0based + 1, end_1based, len(motif), int((end_1based - start_0based)/len(motif)),
                        f"{chrom}-{start_0based}-{end_1based}-{motif}"
                    ]
                else:
                    raise ValueError("Neither gangstr nor hipstr set to True")

                f.write("\t".join(map(str, output_fields)) + "\n")

    print(f"Wrote {len(batches):,d} GangSTR repeat spec bed files to {output_path_prefix}*.bed")


def write_bed_files(locus_set, output_path):
    with open(os.path.expanduser(output_path), "wt") as f:
        for chrom, start_0based, end_1based, motif in sorted(locus_set):
            f.write("\t".join(map(str, [chrom, start_0based, end_1based, motif, len(motif)])) + "\n")

    os.system(f"bgzip -f {output_path}")
    os.system(f"tabix -f {output_path}.gz")
    print(f"Wrote {len(locus_set):,d} loci to {output_path}.gz")


def main():
    args = parse_args()

    if re.search(".bed(.gz)?$", args.truth_set_variants_tsv_or_bed_path):
        truth_set_df = pd.read_table(args.truth_set_variants_tsv_or_bed_path,
                                     names=["Chrom", "Start0Based", "End1Based", "Motif", "_"],
                                     dtype={"Chrom": str, "Start0Based": int, "End1Based": int, "Motif": str})

        truth_set_df.loc[:, "Start1Based"] = truth_set_df.Start0Based + 1
        truth_set_df.loc[:, "IsFoundInReference"] = "Yes"
    else:
        # Read in STR truth set
        truth_set_df = pd.read_table(args.truth_set_variants_tsv_or_bed_path)

    print(f"Parsed {len(truth_set_df):,d} records from {args.truth_set_variants_tsv_or_bed_path}")

    # Must be found in reference for EH and GangSTR to work
    length_before = len(truth_set_df)
    truth_set_df = truth_set_df[truth_set_df["IsFoundInReference"] == "Yes"]
    if length_before > len(truth_set_df):
        print(f"Discarded {length_before - len(truth_set_df):,d} loci without matching repeats in the reference")

    # Generate positive (ie. variant) and negative (non-variant) loci for the variant catalogs
    positive_loci, positive_loci_counters, truth_set_loci_interval_trees = generate_set_of_positive_loci(truth_set_df)
    print(f"Generated {len(positive_loci):,d} positive loci")

    locus_sets = [("positive", positive_loci)]
    if not args.skip_negative_loci:
        all_pure_repeats_in_syndip_high_confidence_regions = load_all_pure_repeats_in_syndip_confidence_regions(
            args.all_repeats_bed, args.high_confidence_regions_bed)
        print(f"Loaded {len(all_pure_repeats_in_syndip_high_confidence_regions):,d} STRs in hg38 within the SynDip high-confidence regions")

        negative_loci = generate_set_of_negative_loci(
            all_pure_repeats_in_syndip_high_confidence_regions, truth_set_loci_interval_trees, positive_loci_counters)
        print(f"Generated {len(negative_loci):,d} negative loci")
        locus_sets.append(("negative", negative_loci))

    # Generate variant catalogs
    subdirs_to_create = []
    if not args.skip_eh_catalog:       subdirs_to_create.append("expansion_hunter")
    if not args.skip_gangstr_catalog:  subdirs_to_create.append("gangstr")
    if not args.skip_hipstr_catalog:   subdirs_to_create.append("hipstr")

    output_dir = args.output_dir
    for label, locus_set in locus_sets:
        for subdir in subdirs_to_create:
            subdir_path = os.path.join(output_dir, subdir)
            if not os.path.isdir(subdir_path):
                print(f"Creating directory {subdir_path}")
                os.makedirs(subdir_path)

        if not args.skip_eh_catalog:
            write_expansion_hunter_variant_catalogs(locus_set,
                os.path.join(output_dir, f"expansion_hunter/{label}_loci.EHv5"),
                loci_per_run=args.expansion_hunter_loci_per_run)
        if not args.skip_gangstr_catalog:
            write_gangstr_or_hipstr_repeat_specs(locus_set,
                 os.path.join(output_dir, f"gangstr/{label}_loci.GangSTR"),
                 gangstr=True, loci_per_run=args.gangstr_loci_per_run)
        if not args.skip_hipstr_catalog:
            write_gangstr_or_hipstr_repeat_specs(locus_set,
                 os.path.join(output_dir, f"hipstr/{label}_loci.HipSTR"),
                 hipstr=True, loci_per_run=args.gangstr_loci_per_run)
        write_bed_files(locus_set, os.path.join(output_dir, f"{label}_loci.bed"))

    # Make sure positive regions and negative regions don't overlap.
    positive_loci_bedtool = pybedtools.BedTool(
        os.path.expanduser(os.path.join(output_dir, f"positive_loci.bed.gz")))
    if not args.skip_negative_loci:
        overlap_count = positive_loci_bedtool.intersect(
            os.path.expanduser(os.path.join(output_dir, f"negative_loci.bed.gz"))).count()

        if overlap_count > 0:
            raise ValueError(f"ERROR: {overlap_count} loci overlap between positive set and negative set")

    print("Done")


if __name__ == "__main__":
    main()
