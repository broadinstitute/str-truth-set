"""
Create ExpansionHunter and GangSTR variant catalogs for all loci in the STR truth set as well as an equal
number of loci from other places in the genome to use as true-negatives.
"""

import argparse
import collections
import gzip
import json
import os
import pandas as pd
import pybedtools
import pyfaidx
import random
import re

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif

"""This threshold defines how far away an STR locus must be from any indels in the syndip truth set before it is 
considered a negative (ie. non-variant locus).
"""

MIN_DISTANCE_TO_INDELS_AROUND_NEGATIVE_LOCI = 100  # base pairs


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("-R", "--ref-fasta", help="Reference fasta path")
    p.add_argument("--expansion-hunter-loci-per-run", type=int, default=1000, help="ExpansionHunter batch size. "
                   "The set of all STR loci in the truth set will be split into variant catalogs of this size.")
    p.add_argument("--gangstr-loci-per-run", type=int, default=100000, help="GangSTR batch size. "
                   "The set of all STR loci in the truth set will be split into repeat spec files of this size.")
    p.add_argument("--straglr-loci-per-run", type=int, default=10000, help="Straglr batch size. "
                   "The set of all STR loci in the truth set will be split into bed files of this size.")
    p.add_argument("--truth-set-bed", default="./STR_truth_set.v1.variants.bed.gz")
    p.add_argument("--skip-negative-loci", action="store_true", help="Also create a catalog generating negative loci")
    p.add_argument("--high-confidence-regions-bed", default="./ref/full.38.bed.gz",
                   help="Path of the SynDip high-confidence regions .bed file")
    p.add_argument("--all-indels-vcf", #default="./ref/full.38.INDELs.vcf.gz",
                   help="Path of the original SynDip vcf filtered to INDEL variants")
    p.add_argument("--all-hg38-repeats-bed", default="./ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_9bp.bed.gz",
                   help="Path of bed file containing all repeats in the reference genome generated using a tool like "
                        "TandemRepeatFinder")
    p.add_argument("--only", choices=["eh", "gangstr", "hipstr", "popstr", "trgt", "straglr", "longtr"], action="append",
                   help="Only generate catalogs for the specified tool(s)")
    p.add_argument("--skip-eh", action="store_true", help="Skip generating an ExpansionHunter catalog")
    p.add_argument("--skip-gangstr", action="store_true", help="Skip generating a GangSTR catalog")
    p.add_argument("--skip-hipstr", action="store_true", help="Skip generating a HipSTR catalog")
    p.add_argument("--skip-popstr", action="store_true", help="Skip generating a popSTR catalog")
    p.add_argument("--skip-trgt", action="store_true", help="Skip generating a TRGT catalog")
    p.add_argument("--skip-straglr", action="store_true", help="Skip generating a Straglr catalog")
    p.add_argument("--skip-longtr", action="store_true", help="Skip generating a LongTR catalog")
    p.add_argument("--output-dir", default="./tool_comparison/variant_catalogs/", help="Directory where to write output files")
    p.add_argument("--output-filename-prefix", help="Prefix of output files")
    p.add_argument("truth_set_variants_tsv_or_bed_path", nargs="?", default="STR_truth_set.v1.variants.tsv",
                   help="Path of the STR truth set .variants.tsv or of an arbitrary bed file")
    args = p.parse_args()

    print("Args:")
    for i, (label, path) in enumerate([
        ("truth set variants", args.truth_set_variants_tsv_or_bed_path),
    ] + ([] if args.skip_negative_loci else [
        ("--high-confidence-regions-bed", args.high_confidence_regions_bed),
        ("--all-indels-vcf", args.all_indels_vcf),
        ("--truth-set-bed", args.truth_set_bed),
        ("--all-hg38-repeats-bed", args.all_hg38_repeats_bed),
    ])):
        if path and not os.path.isfile(path):
            p.error(f"{path} not found")
        print(f"{label:>30s}: {path}")

    if not args.skip_popstr and (not args.only or "popstr" in args.only):
        if args.ref_fasta is None:
            p.error("--ref-fasta must be set to write popSTR catalogs")
        if not os.path.isfile(args.ref_fasta):
            p.error(f"{args.ref_fasta} not found")

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
        high_confidence_regions_bed_path,
        all_indels_vcf_path,
        truth_set_loci_bed_path,
        positive_loci_counters):
    """Generate negative loci by starting with all pure repeats in hg38, filtering them to those that are within
    high-confidence regions, and then excluding all that are near InDels (in case these are STRs)
    or are in the STR truth set. Finally, downsample this set of negative loci to a smaller set so that the number of
    negative loci for each motif size is approximately equal to the number of positive loci for that motif size.
    """

    pure_repeats_bedtool = pybedtools.BedTool(os.path.expanduser(all_hg38_repeats_bed_path))

    negative_loci_bedtool = pure_repeats_bedtool.intersect(
        high_confidence_regions_bed_path, f=1, wa=True, u=True)

    if all_indels_vcf_path:
        negative_loci_bedtool = negative_loci_bedtool.window(
            all_indels_vcf_path, w=MIN_DISTANCE_TO_INDELS_AROUND_NEGATIVE_LOCI, v=True)

    negative_loci_bedtool = negative_loci_bedtool.intersect(
        truth_set_loci_bed_path, v=True)

    negative_loci = []
    for bed_row in negative_loci_bedtool:
        negative_loci.append(bed_row[:4])
    random.seed(0)

    print(f"Got {len(negative_loci):,d} total negative loci")

    # shuffle records to avoid selecting negative loci mostly from chr1
    random.shuffle(negative_loci)

    # compute negative loci
    print(f"Downsampling this to match the number of positive loci for each motif size")
    downsampled_negative_loci = []
    negative_loci_counters = collections.defaultdict(int)
    enough_negative_loci = collections.defaultdict(bool)
    for chrom, start_0based, end, motif in negative_loci:
        if "N" in motif.upper():
            continue
        key = compute_row_key(motif)
        if enough_negative_loci[key]:
            continue

        downsampled_negative_loci.append((chrom, int(start_0based), int(end), motif))
        negative_loci_counters[key] += 1
        if negative_loci_counters[key] >= positive_loci_counters[key]:
            enough_negative_loci[key] = True

            # check if finished
            if all(enough_negative_loci.values()):
                break

    num_duplicate_loci = len(downsampled_negative_loci) - len(set(downsampled_negative_loci))
    if num_duplicate_loci > 0:
        raise ValueError(f"{num_duplicate_loci} negative loci are duplicates")

    return downsampled_negative_loci


def write_expansion_hunter_variant_catalogs(locus_set, output_path_prefix, loci_per_run=None):
    variant_catalog = []

    for unmodified_chrom, start_0based, end_1based, motif in sorted(
        # sort by canonical motif to improve cache hit rate for the optimized version of ExpansionHunter
        locus_set, key=lambda x: compute_canonical_motif(x[3], include_reverse_complement=True)
    ):
        chrom = unmodified_chrom.replace("chr", "")
        variant_catalog.append({
            "LocusId": f"{chrom}-{start_0based}-{end_1based}-{motif}",
            "ReferenceRegion": f"{unmodified_chrom}:{start_0based}-{end_1based}",
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


def compute_repeat_purity(ref_sequence, pure_sequence):
    if len(ref_sequence) != len(pure_sequence):
        raise ValueError(f"ref_sequence and pure_sequence must be the same length: "
                         f"{len(ref_sequence)} != {len(pure_sequence)}: {ref_sequence} vs {pure_sequence}")

    match_counter = sum(1 for a, b in zip(ref_sequence, pure_sequence) if a.upper() == b.upper())
    return match_counter / len(ref_sequence)


def write_popstr_catalogs(locus_set, fasta_obj, output_path_prefix):
    """Example row from a popSTR markerInfo file:
$1                     chrom : chr22
$2           startCoordinate : 10517060  (1-based)
$3             endCoordinate : 10517070
$4               repeatMotif : ATGAG
$5         numOfRepeatsInRef : 2.2
$6   1000refBasesBeforeStart : AAAAGTCCATCATCAAATGAACAGATGAAGAAAATATGGTATATGTGTGTGTGTGTATATATATATATGTATGGTATATATATGTATGGTGTATATATATATATATGTATGATATATATAGTATGGCATATATATGTATGGTGTATATATATGGTATATTTATGTGTATATATATGTATATATGTATATATATGTATATATACATACACACACAGAATGGAATATTAGTCAGCCTTCAAAAGGAAAATTCTGTCGTATTTCAACATGTATCAAGCTTAAGGATATTGTGCTAAGTGAAATAAGCCAGACACAAAGACAAATATATCGTGATTCCATTTATATGATGTATTGAAAGTAGCCAAACACATGGAAACACAAGATAAAATGGTAGTGGTCAGGGCCTGGAGGAAACAGGAAATCTGGAGTTGCTGTTCACCAGGTGTAGAGTTTCAGTCAAGCAAGATAAAAACATTCTAGATTTCTGCTGTACAACAATGTGTATATCATTAACAAAATGTTCTGCAAACTTAAAATTTTGTTAAGATGGTAGATTTTTTTGTTATGTGTTTTTTAATTACAAAAATTTCTGTCTGTATTTAGTTTACATTTTAGTAAGAAAAGACAAATAACCTAATAAATGGTAATGATAAATGCTGTAAGGATATCTAGAGCAATAAAAGAAGTTATGGATATGGGAGTATAATTTTAGATAGTGAAGATTTCTGTATTCAAATGCCACATGCAAAAAGGACTAAGGGGAATGAGGAGATGAGTCATATGGAATGCCCAAGACACAGAAAAAGGCAGGCAGAGAAAATAACAAGGATAAAGACACTGAAGTAAAATCATGCTTTCTATGTTTCAAAACAGCAGCAAGTACACTACTGTGATAGAGCAGGAGTGACCAATGGGAGGCTGGAGATTTTGTCAGAGATATTGTCAAGGCTCAGGATCGTACAGGGACTTGTAAGCCTGGAAAGCACTTTGAATTTTATTCAGA
$7      1000refBasesAfterEnd : AGCCATTGAAAAGTTTTTAAGCAGATGAGTAAAATAATCCACCTTGTATTTTAAGAGGAGCATTCTACCTTCTCTGTGGAATAGAGAGGTGGAAGGGAAAGCTTGAAGCAGAGAGAGCAGTGAAGAGTGTGCTGTAATATTCTTATGTGAGAAATAGTGGAGGGAATGAGAGGTGGTCAGCCTTAAACTGCCATTTGCTCTCTGTATCAGGGCTCAGGGACTTTCAGACTCTCCAGGGATTCCTTACAGTTTTTACATTTGTTTCTCAGTTGCAAACATAATCTCTTCACTCTATGAGAACTCTAGATCTGAATCCTTGTTATGAGTCAGGAGTCCTCTCTAGTTTACTCTCTTGTCAGTAACTAGACTTGAAATGTTTTGATATTAATTGATATAGTAATAAGATTGGTTAGAGAAATAGCAAAGAGAGAGCATCCCCATCCTATGACCATATCAGCACCAGAAGAGAAAAACACATCTACACAGTTTTTCCCTTGGCATAGGTCCTGGTATTCTGTTAGGGAACAGGTTATGAAGGAAATACAAAGGGTCTTTGTACTTACTTTCAGAATGTATTTTCTTTAACATGAAAAGAATCCAAGGCCTTTTTGCTTCTAATTGCTTTTTGTGTATCTACTACCATCCCTTGCTAAATTATTGATAGGTTTCCTCAAATCTCGGCATGATGTCCTACATTCTAAATTTTCAATAGCTGAAAATTTCACCTTTTCAGTGCCTTCAAGTTTATCTCAGTAAAAAGTTGAGAAAGACTGTAATAGAGTTATTTAATCAGATTTTTTTCATCTACCATAATTTTTGAATAAGGAAAAACAGCAATACTTTTTCTCCTTACTTGGCAAGTAATTTTCATAGAGAGGAAAAAAACAATCAAAACAGGTACAAAATGTAACAAAACCAAAGGACCATGTGAGGTGAAATTTAAAATGAGAAAAATGTCCACAGTACTTTGGGCAATGCAACTCCTGAGAAATAGTAACTC
$8          repeatSeqFromRef : ATGAGATGAGA
$9              minFlankLeft : 4
$10            minFlankRight : 4
$11             repeatPurity : 1.00
$12          defaultSlippage : 0.017275
$13           defaultStutter : 0.946913
$14         fractionAinMotif : 0.4
$15         fractionCinMotif : 0
$16         fractionGinMotif : 0.4
$17         fractionTinMotif : 0.2
    """
    f = None
    previous_chrom = None
    for chrom, start_0based, end_1based, motif in sorted(locus_set):
        if len(motif) > 6:
            continue

        if not chrom.startswith("chr"):
            chrom = f"chr{chrom}"

        if chrom != previous_chrom:
            previous_chrom = chrom
            output_path = f"{output_path_prefix}.{chrom}.markerInfo.gz"
            if f is not None:
                f.close()
            f = gzip.open(output_path, "wt")
            print(f"Writing to {output_path}")

        repeat_seq_from_ref = fasta_obj[chrom][start_0based : end_1based]
        num_repeats_in_ref = len(repeat_seq_from_ref) / len(motif)
        pure_repeat_seq = motif * int(num_repeats_in_ref + 1)
        f.write(" ".join(map(str, [
            chrom,                  # chrom
            str(start_0based + 1),  # startCoordinate
            str(end_1based),        # endCoordinate
            motif,                  # repeatMotif
            "%0.1f" % num_repeats_in_ref,                            # numOfRepeatsInRef
            fasta_obj[chrom][start_0based - 1000 : start_0based],   # 1000refBasesBeforeStart
            fasta_obj[chrom][end_1based : end_1based + 1000],       # 1000refBasesAfterEnd
            repeat_seq_from_ref,                                    # repeatSeqFromRef
            4,     # minFlankLeft
            4,     # minFlankRight
            compute_repeat_purity(repeat_seq_from_ref, pure_repeat_seq[:len(repeat_seq_from_ref)]),  # repeatPurity
            0.02,  # defaultSlippage
            0.95,  # defaultStutter
            "%0.1f" % (motif.upper().count("A") / len(motif)),  # fractionAinMotif
            "%0.1f" % (motif.upper().count("C") / len(motif)),  # fractionCinMotif
            "%0.1f" % (motif.upper().count("G") / len(motif)),  # fractionGinMotif
            "%0.1f" % (motif.upper().count("T") / len(motif)),  # fractionTinMotif
        ])) + "\n")

    if f is not None:
        f.close()


def write_gangstr_hipstr_or_longtr_repeat_specs(locus_set, output_path_prefix, tool="gangstr", loci_per_run=None):
    locus_list = list(locus_set)

    if tool not in ("gangstr", "hipstr", "longtr"):
        raise ValueError(f"Invalid tool arg: '{tool}'. Must be 'gangstr', 'hipstr', or 'longtr'")

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
                if end_1based - start_0based <= 1:
                    # This can happen with homopolymers where there is only 1 repeat in the reference, so it passes
                    # the IsFoundInReference filter, but then errors out for HipSTR, and LongTR
                    continue
                chrom = chrom
                if tool == "gangstr":
                    output_fields = [chrom, start_0based + 1, end_1based, len(motif), motif]
                elif tool == "longtr":
                    output_fields = [
                        chrom, start_0based + 1, end_1based, len(motif), int((end_1based - start_0based)/len(motif)),
                        f"{chrom}-{start_0based}-{end_1based}-{motif}"
                    ]
                elif tool == "hipstr":
                    if len(motif) > 9:
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

    print(f"Wrote {len(batches):,d} {tool} repeat spec bed files to {output_path_prefix}*.bed")


def write_trgt_catalog(locus_set, output_path):
    with open(os.path.expanduser(output_path), "wt") as f:
        # example: chr1	1224281	1224291	ID=chr1_1224281_1224291;MOTIFS=TTTTA;STRUC=(TTTTA)n
        for chrom, start_0based, end_1based, motif in sorted(locus_set):
            column4 = f"ID={chrom}_{start_0based}_{end_1based};MOTIFS={motif};STRUC=({motif})n"
            f.write("\t".join(map(str, [chrom, start_0based, end_1based, column4])) + "\n")

    print(f"Wrote {len(locus_set):,d} loci to {output_path}")


def write_straglr_catalog(locus_set, output_path_prefix, loci_per_run=None):
    if loci_per_run is None:
        batches = [locus_set]
    else:
        locus_list = list(locus_set)
        batches = [
            locus_list[i:i+loci_per_run] for i in range(0, len(locus_list), loci_per_run)
        ]

    for batch_i, current_locus_list in enumerate(batches):
        with open(f"{output_path_prefix}.{batch_i+1:03d}_of_{len(batches):03d}.bed", "wt") as f:
            for chrom, start_0based, end_1based, motif in sorted(current_locus_list):
                f.write("\t".join(map(str, [chrom, start_0based, end_1based, motif])) + "\n")

    print(f"Wrote {len(batches):,d} Straglr catalogs to {output_path_prefix}*.bed")


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
        print("Converting bed file to tsv", args.truth_set_variants_tsv_or_bed_path)
        truth_set_df = pd.read_table(args.truth_set_variants_tsv_or_bed_path, index_col=False,
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
            args.high_confidence_regions_bed,
            args.all_indels_vcf,
            args.truth_set_bed,
            positive_loci_counters)

        print(f"Generated {len(negative_loci):,d} negative loci")
        locus_sets.append(("negative", negative_loci))

    # Generate variant catalogs
    subdirs_to_create = []
    if not args.skip_eh and (not args.only or "eh" in args.only): subdirs_to_create.append("expansion_hunter")
    if not args.skip_gangstr and (not args.only or "gangstr" in args.only): subdirs_to_create.append("gangstr")
    if not args.skip_hipstr and (not args.only or "hipstr" in args.only):   subdirs_to_create.append("hipstr")
    if not args.skip_popstr and (not args.only or "popstr" in args.only):   subdirs_to_create.append("popstr")
    if not args.skip_trgt and (not args.only or "trgt" in args.only):   subdirs_to_create.append("trgt")
    if not args.skip_straglr and (not args.only or "straglr" in args.only):   subdirs_to_create.append("straglr")
    if not args.skip_longtr and (not args.only or "longtr" in args.only):   subdirs_to_create.append("longtr")

    fasta_obj = None
    output_dir = args.output_dir
    output_filename_prefix = f"{args.output_filename_prefix}." if args.output_filename_prefix else ""

    for label, locus_set in locus_sets:
        for subdir in subdirs_to_create:
            subdir_path = os.path.join(output_dir, subdir)
            if not os.path.isdir(subdir_path):
                print(f"Creating directory {subdir_path}")
                os.makedirs(subdir_path)

        if not args.skip_eh and (not args.only or "eh" in args.only):
            write_expansion_hunter_variant_catalogs(locus_set,
                os.path.join(output_dir, f"expansion_hunter/{output_filename_prefix}{label}_loci.EHv5"),
                loci_per_run=args.expansion_hunter_loci_per_run)

            # also output a single catalog with all loci
            if args.expansion_hunter_loci_per_run < len(locus_set):
                write_expansion_hunter_variant_catalogs(locus_set,
                    os.path.join(output_dir, f"expansion_hunter/{output_filename_prefix}{label}_loci.EHv5"),
                    loci_per_run=10**9)

        if not args.skip_gangstr and (not args.only or "gangstr" in args.only):
            write_gangstr_hipstr_or_longtr_repeat_specs(locus_set,
                 os.path.join(output_dir, f"gangstr/{output_filename_prefix}{label}_loci.GangSTR"),
                 tool="gangstr", loci_per_run=args.gangstr_loci_per_run)

        if not args.skip_hipstr and (not args.only or "hipstr" in args.only):
            write_gangstr_hipstr_or_longtr_repeat_specs(locus_set,
                 os.path.join(output_dir, f"hipstr/{output_filename_prefix}{label}_loci.HipSTR"),
                 tool="hipstr", loci_per_run=args.gangstr_loci_per_run)

        if not args.skip_popstr and (not args.only or "popstr" in args.only):
            if fasta_obj is None:
                fasta_obj = pyfaidx.Fasta(args.ref_fasta, one_based_attributes=False, as_raw=True, sequence_always_upper=True)

            write_popstr_catalogs(locus_set, fasta_obj,
                 os.path.join(output_dir, f"popstr/{output_filename_prefix}{label}_loci.popSTR"))

        if not args.skip_trgt and (not args.only or "trgt" in args.only):
            write_trgt_catalog(locus_set, os.path.join(output_dir, f"trgt/{output_filename_prefix}{label}_loci.TRGT_repeat_catalog.bed"))

        if not args.skip_straglr and (not args.only or "straglr" in args.only):
            write_straglr_catalog(locus_set, os.path.join(output_dir, f"straglr/{output_filename_prefix}{label}_loci.straglr_catalog"),
                                  loci_per_run=args.straglr_loci_per_run)

        if not args.skip_longtr and (not args.only or "longtr" in args.only):
            write_gangstr_hipstr_or_longtr_repeat_specs(locus_set,
                os.path.join(output_dir, f"longtr/{output_filename_prefix}{label}_loci.LongTR"),
                tool="longtr", loci_per_run=None)

        write_bed_files(locus_set, os.path.join(output_dir, f"{output_filename_prefix}{label}_loci.bed"))

    # Make sure positive regions and negative regions don't overlap.
    positive_loci_bedtool = pybedtools.BedTool(
        os.path.expanduser(os.path.join(output_dir, f"{output_filename_prefix}positive_loci.bed.gz")))
    if not args.skip_negative_loci:
        overlap_count = positive_loci_bedtool.intersect(
            os.path.expanduser(os.path.join(output_dir, f"{output_filename_prefix}negative_loci.bed.gz"))).count()

        if overlap_count > 0:
            raise ValueError(f"ERROR: {overlap_count} loci overlap between positive set and negative set")

    print("Done")


if __name__ == "__main__":
    main()
