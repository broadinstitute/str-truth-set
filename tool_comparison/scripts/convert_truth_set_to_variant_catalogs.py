"""
Create ExpansionHunter and GangSTR variant catalogs for all loci in the STR truth set as well as an equal
number of loci from other places in the genome to use as true-negatives.
"""

import argparse
import collections
import json
import os
import pandas as pd
import pybedtools
import random
import re

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


"""This threshold defines how far away an STR locus must be from any indels in the syndip truth set before it is 
considered a negative (ie. non-variant locus).
"""

MIN_DISTANCE_TO_INDELS_AROUND_NEGATIVE_LOCI = 100  # base pairs


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--expansion-hunter-loci-per-run", type=int, default=500, help="ExpansionHunter batch size. "
                   "The set of all STR loci in the truth set will be split into variant catalogs of this size.")
    p.add_argument("--gangstr-loci-per-run", type=int, default=10000, help="GangSTR batch size. "
                   "The set of all STR loci in the truth set will be split into repeat spec files of this size.")
    p.add_argument("--syndip-high-confidence-regions-bed", default="./ref/full.38.bed.gz",
                    help="Path of the SynDip high-confidence regions .bed file")
    p.add_argument("--syndip-indels-vcf", default="./ref/full.38.INDELs.vcf.gz",
                   help="Path of the original SynDip vcf")
    p.add_argument("--truth-set-bed", default="./STR_truth_set.v1.variants.bed.gz")
    p.add_argument("--all-hg38-repeats-bed", default="./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz",
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

    for i, path in enumerate([
        args.truth_set_variants_tsv_or_bed_path,
        args.syndip_high_confidence_regions_bed,
        args.syndip_indels_vcf,
        args.truth_set_bed,
        args.all_hg38_repeats_bed,
    ]):
        if not os.path.isfile(path):
            p.error(f"{path} not found")
        print(f"Input {i+1}: {path}")
    return args


def compute_row_key(motif):
    if len(motif) <= 4:
        # stratify loci only by their canonical motif
        canonical_motif = compute_canonical_motif(motif, include_reverse_complement=True)
        key = (len(canonical_motif), canonical_motif)
    elif len(motif) <= 24:
        # above 4 repeats, stratify loci only by motif size and not the motif itself
        key = (len(motif), None)
    else:
        # group motifs of size 25bp or larger into one bin
        key = (25, None)

    return key


def generate_set_of_positive_loci(truth_set_df):
    """Convert the truth set table into a set of positive loci, as well as a dictionary of counters stratified by motif,
    and interval trees for efficiently computing interval overlaps with loci in the truth set.
    """

    positive_loci = set()  # loci with STR variants in CHM1-CHM13
    positive_loci_counters = collections.defaultdict(int)
    for _, row in truth_set_df.iterrows():
        motif = row.Motif
        chrom = row.Chrom
        start_0based = int(row.Start1Based) - 1
        end = int(row.End1Based)
        if (end - start_0based) % len(motif) != 0:
            raise ValueError(f"{row.to_dict()} reference interval size is not a multiple of the motif")

        key = compute_row_key(motif)
        positive_loci.add((chrom, start_0based, end, motif))
        positive_loci_counters[key] += 1

    return positive_loci, positive_loci_counters


def generate_set_of_negative_loci(
        all_hg38_repeats_bed_path,
        syndip_high_confidence_regions_bed_path,
        syndip_indels_vcf_path,
        truth_set_loci_bed_path,
        positive_loci_counters):

    pure_repeats_bedtool = pybedtools.BedTool(os.path.expanduser(all_hg38_repeats_bed_path))
    negative_loci_bedtool = pure_repeats_bedtool.intersect(syndip_high_confidence_regions_bed_path, f=1, wa=True, u=True)
    negative_loci_bedtool = negative_loci_bedtool.window(syndip_indels_vcf_path, w=MIN_DISTANCE_TO_INDELS_AROUND_NEGATIVE_LOCI, v=True)
    negative_loci_bedtool = negative_loci_bedtool.intersect(truth_set_loci_bed_path, v=True)

    negative_loci = []
    for bed_row in negative_loci_bedtool:
        negative_loci.append(bed_row[:4])
    random.seed(0)

    # shuffle records to avoid selecting negative loci mostly from chr1
    random.shuffle(negative_loci)

    # compute negative loci
    negative_loci = set()
    negative_loci_counters = collections.defaultdict(int)
    enough_negative_loci = collections.defaultdict(bool)
    for chrom, start_0based, end, motif in negative_loci:
        key = compute_row_key(motif)
        if enough_negative_loci[key]:
            continue

        negative_loci.add((chrom, int(start_0based), int(end), motif))
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
                chrom = chrom
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
        truth_set_df.loc[:, "IsFoundInReference"] = True
    else:
        # Read in STR truth set
        truth_set_df = pd.read_table(args.truth_set_variants_tsv_or_bed_path)

    print(f"Parsed {len(truth_set_df):,d} records from {args.truth_set_variants_tsv_or_bed_path}")

    # Must be found in reference for EH and GangSTR to work
    length_before = len(truth_set_df)
    truth_set_df = truth_set_df[truth_set_df["IsFoundInReference"]]
    if length_before > len(truth_set_df):
        print(f"Discarded {length_before - len(truth_set_df):,d} loci without matching repeats in the reference")

    # Generate positive (ie. variant) loci for the variant catalogs
    positive_loci, positive_loci_counters = generate_set_of_positive_loci(truth_set_df)
    print(f"Generated {len(positive_loci):,d} positive loci")

    locus_sets = [("positive", positive_loci)]

    if not args.skip_negative_loci:
        # generate negative (non-variant) loci for the variant catalogs
        negative_loci = generate_set_of_negative_loci(
            args.all_hg38_repeats_bed,
            args.syndip_high_confidence_regions_bed,
            args.syndip_indels_vcf,
            args.truth_set_bed,
            positive_loci_counters)

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
