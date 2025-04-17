"""This script converts a ConSTRain output VCF to the ExpansionHunter output .json format
which I use as a common input format in downstream scripts.
"""


"""
ExpansionHunter output format:

  "LocusResults": {
        "chr12-57610122-57610131-GCA": {
          "AlleleCount": 2,
          "Coverage": 50.469442942130875,
          "FragmentLength": 433,
          "LocusId": "chr12-57610122-57610131-GCA",
          "ReadLength": 151,
          "Variants": {
            "chr12-57610122-57610131-GCA": {
              "CountsOfFlankingReads": "(1, 1), (2, 4)",
              "CountsOfInrepeatReads": "()",
              "CountsOfSpanningReads": "(2, 1), (3, 48), (6, 1)",
              "Genotype": "3/3",
              "GenotypeConfidenceInterval": "3-3/3-3",
              "ReferenceRegion": "chr12:57610122-57610131",
              "RepeatUnit": "GCA",
              "VariantId": "chr12-57610122-57610131-GCA",
              "VariantType": "Repeat"
            }
          }
        },

  "SampleParameters": {
        "SampleId": "NA19239",
        "Sex": "Female"
  }
"""


import argparse
import gzip
import simplejson as json
import re
from tqdm import tqdm


def main():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--discard-hom-ref", action="store_true", help="Discard hom-ref calls")
    p.add_argument("--output-REF-ALT-fields", action="store_true", help="Output the VCF REF and ALT fields in the output")
    p.add_argument("--verbose", action="store_true", help="Print verbose output")
    p.add_argument("--show-progress-bar", action="store_true", help="Show a progress bar")
    p.add_argument("--sample-id",
                   help="If not specified, the sample id will be parsed from the last column of the vcf header.")
    p.add_argument("vcf_path", help="TRGT vcf path")
    args = p.parse_args()

    print(f"Processing {args.vcf_path}")
    locus_results = process_constrain_vcf(
        args.vcf_path,
        sample_id=args.sample_id,
        discard_hom_ref=args.discard_hom_ref,
        output_REF_ALT_fields=args.output_REF_ALT_fields,
        verbose=args.verbose,
        show_progress_bar=args.show_progress_bar,
    )

    output_json_path = re.sub(".vcf(.gz)?$", "", args.vcf_path) + ".json"
    print(f"Writing {len(locus_results['LocusResults']):,d} loci to {output_json_path}")
    with open(output_json_path, "wt") as f:
        json.dump(locus_results, f, indent=3, ignore_nan=True)


def process_constrain_vcf(
    vcf_path, sample_id=None, discard_hom_ref=True, output_REF_ALT_fields=False, verbose=False, show_progress_bar=False):

    locus_results = {
        "LocusResults": {},
        "SampleParameters": {
            "SampleId": sample_id,
            "Sex": None,
        },
    }

    fopen = gzip.open if vcf_path.endswith("gz") else open

    with fopen(vcf_path, "rt") as vcf:
        for line in vcf:
            if line.startswith("#CHROM"):
                header_fields = line.strip().split("\t")
                if sample_id is None and len(header_fields) == 10:
                    print(f"Got sample id '{header_fields[9]}' from the VCF header")
                    locus_results["SampleParameters"]["SampleId"] = header_fields[9]
            if not line.startswith("#"):
                break

    with fopen(vcf_path, "rt") as vcf:
        line_counter = 0
        if show_progress_bar:
            vcf = tqdm(vcf, unit=" vcf records", unit_scale=True, unit_divisor=1000)

        for line in vcf:
            if line.startswith("#"):
                continue

            line_counter += 1

            try:
                fields = line.strip().split("\t")
                chrom = fields[0]
                start_1based = int(fields[1])
                info = fields[7]
                if not fields[9] or fields[9] == ".":  # no genotype
                    continue

                info_dict = dict([key_value.split("=") for key_value in info.split(";")])

                genotype_fields = fields[8].split(":")
                genotype_values = fields[9].split(":")
                genotype_dict = dict(zip(genotype_fields, genotype_values))

                if genotype_dict["GT"] == ".":
                    continue
                # GT:FT:CN:DP:FREQS:REPLEN
                if discard_hom_ref and genotype_dict["GT"] == "0/0":
                    continue

                end_1based = int(info_dict["END"])

                motif = info_dict["RU"]
                allele_sizes_bp = [int(num_repeats) * len(motif) for num_repeats in genotype_dict["REPLEN"].split(",")]
                flip_alleles = len(allele_sizes_bp) == 2 and allele_sizes_bp[0] > allele_sizes_bp[1]
                if flip_alleles:
                    for key in "REPLEN", "FREQS":
                        genotype_dict[key] = ",".join(genotype_dict[key].split(",")[::-1])

                locus_id = f"{chrom}-{start_1based - 1}-{end_1based}-{motif}"

                motif_counts = [num_repeats for num_repeats in genotype_dict["REPLEN"].split(",")]
                motif_counts_CIs = [f"{num_repeats}-{num_repeats}" for num_repeats in genotype_dict["REPLEN"].split(",")]

                variant_id = locus_id

                variant_record = {
                    "Genotype": "/".join(motif_counts),
                    "GenotypeConfidenceInterval": "/".join(motif_counts_CIs),
                    "ReferenceRegion": f"{chrom}:{start_1based - 1}-{end_1based}",
                    "RefAlleleSize": int(info_dict["REF"]),
                    "RepeatUnit": motif,
                    "VariantId": variant_id,
                    #"VariantType": "Repeat",
                }

                if output_REF_ALT_fields:
                    variant_record["Ref"] = fields[3]
                    variant_record["Alt"] = fields[4]

                locus_results["LocusResults"][locus_id] = {
                    "AlleleCount": 2,
                    "LocusId": locus_id,
                    "Variants": {
                        variant_id: variant_record,
                    },
                }

            except Exception as e:
                print(f"Error while parsing vcf record #{line_counter}: {e}")
                print(f"    {line}")
                import traceback
                traceback.print_exc() # print stack trace


    return locus_results


if __name__ == "__main__":
    main()
