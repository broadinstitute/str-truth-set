import argparse
import collections
import gzip
from numbers_utils import search
import os
import pandas as pd

from str_analysis.filter_vcf_to_STR_variants import compute_indel_variant_bases


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detailed", help="Include detailed stats", action="store_true")
    parser.add_argument("--high-confidence-vcf", help="Path of VCF", default="step1.high_confidence_regions.vcf.gz")
    parser.add_argument("--alleles-table", help="Path of alleles table", default="step2.STRs_including_homopolymers.alleles.tsv.gz")
    parser.add_argument("--filtered-out-indels-vcf", help="Path of VCF", default="step2.STRs_including_homopolymers.filtered_out_indels.vcf.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",  default="figures_and_tables/table4_indel_stats.html")
    args = parser.parse_args()

    if args.detailed:
        args.output_html = args.output_html.replace(".html", ".detailed.html")

    # check arg files exist
    for arg_name, arg_value in sorted(vars(args).items()):
        if arg_name.endswith("_table") or arg_name.endswith("_vcf"):
            if not os.path.exists(arg_value):
                parser.error(f"File does not exist: {arg_value}")

    # print arg values
    print("Args:")
    for arg_name, arg_value in sorted(vars(args).items()):
        print(f"{arg_name:30s} {arg_value}")

    with open(args.step_A_log, "rt") as f:
        stepA_log_contents = f.read()

    high_confidence_alleles_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)
    high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
    high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))
    high_confidence_indel_alleles_from_step_A_log = high_confidence_INS_alleles + high_confidence_DEL_alleles

    high_confidence_indel_alleles_from_vcf = 0
    high_confidence_indel_alleles_from_vcf_multiallelic = 0
    with gzip.open(args.high_confidence_vcf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            ref_allele = fields[3]
            alt_alleles = [alt for alt in fields[4].split(",") if alt != "*"]
            for alt_allele in alt_alleles:
                if len(alt_allele) == len(ref_allele):
                    continue

                variant_bases = compute_indel_variant_bases(ref_allele, alt_allele)
                if variant_bases is None:
                    continue
                high_confidence_indel_alleles_from_vcf += 1
                if len(alt_alleles) > 1:
                    high_confidence_indel_alleles_from_vcf_multiallelic += 1


    if high_confidence_indel_alleles_from_vcf != high_confidence_indel_alleles_from_step_A_log:
        raise ValueError(f"high_confidence_indel_alleles_from_vcf != high_confidence_indel_alleles_from_step_A_log: "
                         f"{high_confidence_indel_alleles_from_vcf:,d} != {high_confidence_indel_alleles_from_step_A_log:,d}")

    df_TR_alleles = pd.read_table(args.alleles_table)
    df_homopolymers = df_TR_alleles[df_TR_alleles.MotifSize == 1]
    df_TR_alleles  = df_TR_alleles[df_TR_alleles.MotifSize > 1]

    total_indel_alleles_key = "total indel alleles within SynDip high confidence regions"
    allele_counters = collections.defaultdict(int)
    multiallelic_counters = collections.defaultdict(int)

    allele_counters[total_indel_alleles_key] += high_confidence_indel_alleles_from_vcf
    multiallelic_counters[total_indel_alleles_key] += high_confidence_indel_alleles_from_vcf_multiallelic

    allele_counters["TR alleles that passed all filter criteria (pre-T2T validation)"] += len(df_TR_alleles)
    multiallelic_counters["TR alleles that passed all filter criteria (pre-T2T validation)"] += sum(df_TR_alleles.IsMultiallelic)

    #if args.detailed:
    allele_counters[f"TR homopolymer alleles: poly-A or T"] += sum(df_homopolymers.Motif.isin({"A", "T"}))
    multiallelic_counters[f"TR homopolymer alleles: poly-A or T"] += sum(df_homopolymers[df_homopolymers.Motif.isin({"A", "T"})].IsMultiallelic)
    allele_counters[f"TR homopolymer alleles: poly-C or G"] += sum(df_homopolymers.Motif.isin({"C", "G"}))
    multiallelic_counters[f"TR homopolymer alleles: poly-C or G"] += sum(df_homopolymers[df_homopolymers.Motif.isin({"C", "G"})].IsMultiallelic)
    #else:
    #    allele_counters["TR homopolymer alleles"] += len(df_homopolymers[df_homopolymers.Motif.isin({"A", "C", "G", "T"})])
    #    multiallelic_counters[f"TR homopolymer alleles"] += sum(df_homopolymers[df_homopolymers.Motif.isin({"A", "C", "G", "T"})].IsMultiallelic)

    allele_counters["TR alleles that were excluded for other reasons"] += sum(df_homopolymers.Motif == "N")
    multiallelic_counters["TR alleles that were excluded for other reasons"] += sum(df_homopolymers[df_homopolymers.Motif == "N"].IsMultiallelic)

    with gzip.open(args.filtered_out_indels_vcf, "rt") as filtered_out_indels_vcf:
        for line in filtered_out_indels_vcf:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            ref_allele = fields[3]
            alt_alleles = [alt for alt in fields[4].split(",") if alt != "*"]
            filter_values = fields[6].split(";")
            if len(alt_alleles) != len(filter_values):
                raise ValueError(f"len(alt_alleles) != len(filter_values): "
                                 f"{len(alt_alleles)} != {len(filter_values)}: {fields[4]} {fields[6]}")

            for alt_allele, filter_value in zip(alt_alleles, filter_values):
                if filter_value == "SNV/MNV":
                    # skip since we're only interested in indels
                    continue

                if filter_value == "complex multinucleotide insertion + deletion":
                    # skip variants where multiple bases are inserted and deleted at the same time and which we don't consider simple indels
                    continue

                if filter_value == "repeat unit < 2 bp":
                    raise ValueError(f"Unexpected filter_value: {filter_value} given that input files should include homopolymers")

                if filter_value in {"spans < 9 bp",  "is only 2 repeats",  "INDEL without repeats"}:
                    variant_bases = compute_indel_variant_bases(ref_allele, alt_allele)
                    filter_value = (f"indels with size 1-5bp" if len(variant_bases) <= 5 else "indels with size ≥ 6bp") + f" that did not pass TR filter criteria"
                elif "STR allele" in filter_values:
                    filter_value = f"indel allele in a multi-allelic variant where only one of the alleles passed TR filter criteria"
                else:
                    rename_dict = {
                        "repeat unit > 50 bp":      "TR alleles with pure repeats of motif size > 50bp",
                        "ends in partial repeat":   "TR allele that failed filter because the variant sequence ended in a partial repeat<br /> "
                                                    "Example: chr6:71189213 CAGCAGCA > C",
                        "locus overlaps more than one STR variant":             "TR alleles that were excluded for other reasons",
                        "STR alleles with different motifs":                    "TR alleles that were excluded for other reasons",
                        "STR alleles with different interruption patterns":     "TR alleles that were excluded for other reasons",
                        "STR alleles with different coords":                    "TR alleles that were excluded for other reasons",
                    }
                    if filter_value in rename_dict:
                        filter_value = rename_dict[filter_value]

                allele_counters[filter_value] += 1
                if len(alt_alleles) > 1:
                    multiallelic_counters[filter_value] += 1

    sum_of_counts = sum(allele_counters[k] for k in allele_counters.keys() if k != total_indel_alleles_key)
    if sum_of_counts != high_confidence_indel_alleles_from_step_A_log:
        raise ValueError(f"sum_of_counts != high_confidence_indel_alleles_from_step_A_log: "
                         f"{sum_of_counts:,d} != {high_confidence_indel_alleles_from_step_A_log:,d}")

    print("Allele counts:")
    for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} ({count / high_confidence_indel_alleles_from_step_A_log:6.1%}) {key:100s}  ({multiallelic_counters[key]/count:5.1%} multi-allelic)")

    header = ["Description", "# of indel<br />alleles", "% of indel<br />alleles", "Fraction<br/>multi-allelic"]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"""<tr>{''.join(f'<th nowrap  style="text-align: center;">{h}</th>' for h in header)}</tr>""")

    for i, (key, count) in enumerate(sorted(allele_counters.items(), key=lambda x: -x[1])):
        table_html.append(f"<tr>"
            f"""<td style="line-height: 1.5">{key}</td>"""
            f"""<td style="text-align: right; line-height: 1.5; vertical-align: top">{count:10,d}</td>"""
            f"""<td style="text-align: right; line-height: 1.5; vertical-align: top ">{100 * count / high_confidence_indel_alleles_from_step_A_log:5.1f}%</td>"""
            f"""<td style="text-align: right; line-height: 1.5; vertical-align: top ">{100 * multiallelic_counters[key] / count:5.1f}%</td>"""
        f"</tr>")
        if i == 0:
            table_html.append(f"<tr><td></td><td></td><td></td></tr>")
    table_html.append("</table>")

    with open("figures_and_tables/table_template.html", "rt") as f, open(args.output_html, "wt") as fo:
        table_html_template = f.read()
        table_html_template = table_html_template % {
            "title": "Indel Categories",
            "table": "\n".join(table_html),
        }
        fo.write(table_html_template)


if __name__ == "__main__":
    main()
