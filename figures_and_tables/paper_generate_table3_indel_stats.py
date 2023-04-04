import argparse
import collections
import gzip
from numbers_utils import search
import os
import pandas as pd

from str_analysis.utils.find_repeat_unit import get_most_common_repeat_unit
from str_analysis.filter_vcf_to_STR_variants import compute_indel_variant_bases

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--detailed", help="Include detailed stats", action="store_true")
    parser.add_argument("--alleles-table", help="Path of alleles table", default="step2.STRs.alleles.tsv.gz")
    parser.add_argument("--filtered-out-indels-vcf", help="Path of VCF", default="step2.STRs.filtered_out_indels.vcf.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",  default="figures_and_tables/table3_indel_stats.html")
    args = parser.parse_args()

    if args.detailed:
        args.output_html = args.output_html.replace(".html", ".detailed.html")

    # check arg files exist
    for arg_name, arg_value in sorted(vars(args).items()):
        if arg_name.endswith("_table") or arg_name.endswith("_vcf"):
            if not os.path.exists(arg_value):
                parser.error(f"File does not exist: {arg_value}")

    # print arg values
    for arg_name, arg_value in sorted(vars(args).items()):
        print(f"{arg_name:30s} {arg_value}")

    with open(args.step_A_log, "rt") as f:
        stepA_log_contents = f.read()

    high_confidence_alleles_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)
    high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
    high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))
    high_confidence_indel_alleles_in_syndip = high_confidence_INS_alleles + high_confidence_DEL_alleles

    #high_confidence_indel_alleles_in_syndip = 285701 + 271124

    df_TR_alleles = pd.read_table(args.alleles_table)

    total_indel_alleles_key = "total indel alleles within SynDip high confidence regions"
    allele_counters = collections.defaultdict(int)
    allele_counters[total_indel_alleles_key] = high_confidence_indel_alleles_in_syndip
    allele_counters["TR truth set: TRs that passed all criteria "] = len(df_TR_alleles)

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
                if filter_value in ("SNV/MNV", "complex multinucleotide insertion + deletion"):
                    continue

                #if filter_value == "ends in partial repeat":
                #    print("\t".join(fields[0:5]))

                if filter_value == "INDEL without repeats":
                    if args.detailed:
                        indel_size = abs(len(alt_allele) - len(ref_allele))
                        if indel_size <= 5:
                            filter_value = "indel without repeats: 1bp ≤ indel size ≤ 5bp"
                        elif indel_size <= 10:
                            filter_value = "indel without repeats: 5bp < indel size ≤ 10bp"
                        elif indel_size <= 20:
                            filter_value = "indel without repeats: 10bp < indel size ≤ 20bp"
                        elif indel_size <= 30:
                            filter_value = "indel without repeats: 20bp < indel size ≤ 30bp"
                        elif indel_size > 30:
                            filter_value = "indel without repeats: 30bp < indel size"
                    else:
                        filter_value = "indel without repeats"

                if filter_value == "repeat unit < 2 bp":
                    variant_bases = compute_indel_variant_bases(ref_allele, alt_allele)
                    if variant_bases is None:
                        raise ValueError(f"variant_bases is None: {ref_allele} {alt_allele}")

                    most_common_repeat_unit, most_common_repeat_unit_count = get_most_common_repeat_unit(variant_bases, 1)
                    if args.detailed:
                        if most_common_repeat_unit in ("A", "T"):
                            filter_value = f"homopolymer: poly-A or T"
                        elif most_common_repeat_unit in ("C", "G"):
                            filter_value = f"homopolymer: poly-C or G"
                        elif most_common_repeat_unit == "N":
                            filter_value = "TR allele that failed filter: other reasons"
                        else:
                            raise ValueError(f"Unexpected most_common_repeat_unit: {most_common_repeat_unit}")
                    else:
                        if most_common_repeat_unit == "N":
                            filter_value = "TR allele that failed filter: other reasons"
                        else:
                            filter_value = "homopolymer"

                rename_dict = {
                    "spans < 9 bp":             "TR allele that failed filter criteria: spans < 9bp",
                    "repeat unit > 50 bp":      "TR allele that failed filter criteria: repeat unit > 50bp",
                    "is only 2 repeats":        "TR allele that failed filter criteria: consists of only 2 repeats",
                    "STR allele":               "TR allele passed filter criteria but is part of a multi-allelic <br />"
                                                "variant where the other allele failed filter criteria",
                    "ends in partial repeat":   "TR allele that failed filter: variant sequence ends in partial repeat <br />"
                                                "(Example:  chr6:71189213 CAGCAGCA > C)",
                    "locus overlaps more than one STR variant":             "TR allele that failed filter: other reasons",
                    "STR alleles with different motifs":                    "TR allele that failed filter: other reasons",
                    "STR alleles with different interruption patterns":     "TR allele that failed filter: other reasons",
                    "STR alleles with different coords":                    "TR allele that failed filter: other reasons",
                }
                if filter_value in rename_dict:
                    filter_value = rename_dict[filter_value]

                allele_counters[filter_value] += 1

    sum_of_counts = sum(allele_counters[k] for k in allele_counters.keys() if k != total_indel_alleles_key)
    if sum_of_counts != high_confidence_indel_alleles_in_syndip:
        raise ValueError(f"sum_of_counts != high_confidence_indel_alleles_in_syndip: {sum_of_counts:,d} != {high_confidence_indel_alleles_in_syndip:,d}")

    print("Allele counts:")
    for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} ({100 * count / high_confidence_indel_alleles_in_syndip:5.1f}%) {key:10s}")

    header = ["Description", "# of indels", "% of indels"]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"<tr>{''.join(f'<th nowrap>{h}</th>' for h in header)}</tr>")

    for i, (key, count) in enumerate(sorted(allele_counters.items(), key=lambda x: -x[1])):
        table_html.append(f"<tr>"
            f"""<td style="line-height: 1.5">{key}</td>"""
            f"""<td style="text-align: right; line-height: 1.5; vertical-align: top">{count:10,d}</td>"""
            f"""<td style="text-align: right; line-height: 1.5; vertical-align: top ">{100 * count / high_confidence_indel_alleles_in_syndip:5.1f}%</td>"""
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
