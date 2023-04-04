import argparse
import collections
import gzip
from numbers_utils import search
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alleles-table", help="Path of alleles table", default="step2.STRs.alleles.tsv.gz")
    parser.add_argument("--filtered-out-indels-vcf", help="Path of VCF", default="step2.STRs.filtered_out_indels.vcf.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",  default="figures_and_tables/table3_indel_stats.html")
    args = parser.parse_args()

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

    high_confidence_indel_alleles_in_syndip = 285701 + 271124

    df_TR_alleles = pd.read_table(args.alleles_table)

    allele_counters = collections.defaultdict(int)
    allele_counters["total indel alleles"] = high_confidence_indel_alleles_in_syndip
    allele_counters["total TR alleles in truth set"] = len(df_TR_alleles)

    with gzip.open(args.filtered_out_indels_vcf, "rt") as filtered_out_indels_vcf:
        for line in filtered_out_indels_vcf:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            alt_alleles = [alt for alt in fields[4].split(",") if alt != "*"]
            filter_values = fields[6].split(";")
            if len(alt_alleles) != len(filter_values):
                print(f"WARNING: len(alt_alleles) != len(filter_values): {len(alt_alleles)} != {len(filter_values)}: {fields[4]} {fields[6]}")
            for filter_value in filter_values:
                if filter_value in ("SNV/MNV", "complex multinucleotide insertion + deletion"):
                    continue
                allele_counters[filter_value] += 1

    allele_counters["sum"] = sum(allele_counters[k] for k in allele_counters.keys() if k != "total indel alleles")

    print("Allele counts:")
    for key, count in sorted(allele_counters.items(), key=lambda x: -x[1]):
        print(f"{count:10,d} ({100 * count / high_confidence_indel_alleles_in_syndip:5.1f}%) {key:10s}")

    header = ["Description", "# of indels", "% of indels"]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"<tr>{''.join(f'<th nowrap>{h}</th>' for h in header)}</tr>")

    for i, (key, count) in enumerate(sorted(allele_counters.items(), key=lambda x: -x[1])):
        table_html.append(f"<tr>"
            f"<td>{key}</td>"
            f"<td>{count:10,d}</td>"
            f"<td>{100 * count / high_confidence_indel_alleles_in_syndip:5.1f}%</td>"
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
