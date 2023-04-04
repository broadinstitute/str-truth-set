import argparse
from numbers_utils import search, format_n, format_np
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants-table", help="Path of variants table",  default="step2.STRs.variants.tsv.gz")
    parser.add_argument("--alleles-table",  help="Path of alleles table",   default="step2.STRs.alleles.tsv.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",    default="figures_and_tables/table1_summary_stats.html")
    args = parser.parse_args()

    # check arg files exist
    for arg_name, arg_value in sorted(vars(args).items()):
        if arg_name.endswith("_table") or arg_name.endswith("_log"):
            if not os.path.exists(arg_value):
                parser.error(f"File does not exist: {arg_value}")

    # print arg values
    for arg_name, arg_value in sorted(vars(args).items()):
        print(f"{arg_name:30s} {arg_value}")


    with open(args.step_A_log, "rt") as f:
        stepA_log_contents = f.read()

    df_variants = pd.read_table(args.variants_table)
    df_alleles = pd.read_table(args.alleles_table)

    total_variants_in_syndip = search(f"step1:input:[ ]+([0-9,]+)[ ]+TOTAL[ ]variants", stepA_log_contents, type=int)
    total_alleles_in_syndip = search(f"step1:input:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)

    high_confidence_variants_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]variants", stepA_log_contents, type=int)
    high_confidence_alleles_in_syndip = search(f"step1:output:[ ]+([0-9,]+)[ ]+TOTAL[ ]alleles", stepA_log_contents, type=int)

    high_confidence_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] SNV variants", stepA_log_contents).replace(",", ""))
    high_confidence_multiallelic_SNV_variants = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_variants_in_syndip:,d}.*[)] multiallelic SNV variants", stepA_log_contents).replace(",", ""))
    high_confidence_indel_variants_in_syndip = high_confidence_variants_in_syndip - high_confidence_SNV_variants - high_confidence_multiallelic_SNV_variants

    high_confidence_INS_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] INS alleles", stepA_log_contents).replace(",", ""))
    high_confidence_DEL_alleles = int(search(f"step1:output:[ ]*([0-9,]+) out of[ ]* {high_confidence_alleles_in_syndip:,d}.*[)] DEL alleles", stepA_log_contents).replace(",", ""))
    high_confidence_indel_alleles_in_syndip = high_confidence_INS_alleles + high_confidence_DEL_alleles

    print(f"{format_n(total_variants_in_syndip)} total variants in SynDip")
    print(f"{format_n(total_alleles_in_syndip)} total alleles in SynDip")
    print("")
    print(f"{format_n(high_confidence_variants_in_syndip)} high-confidence variants in SynDip")
    print(f"{format_n(high_confidence_alleles_in_syndip)} high-confidence alleles in SynDip")

    print(f"{format_np(high_confidence_indel_variants_in_syndip, high_confidence_variants_in_syndip)} high-confidence INDEL variants in SynDip")
    print(f"{format_np(high_confidence_indel_alleles_in_syndip, high_confidence_alleles_in_syndip)} high-confidence INDEL alleles in SynDip")

    header = ["Step", "Description", "# of variant loci", "# of alleles", "% of alleles"]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"<tr>{''.join(f'<th nowrap>{h}</th>' for h in header)}</tr>")

    # row 0
    table_html.append("<tr>")
    table_html.append("""<td><i></i></td>""")
    table_html.append("""<td>Start with all variants in SynDip high-confidence regions</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{high_confidence_variants_in_syndip:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{high_confidence_alleles_in_syndip:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">100.0%</td>""")
    table_html.append("</tr>")

    # row 1
    table_html.append("<tr>")
    table_html.append("""<td><i>1</i></td>""")
    table_html.append("""<td>Keep only insertions and deletions</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{high_confidence_indel_variants_in_syndip:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{high_confidence_indel_alleles_in_syndip:,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{high_confidence_indel_alleles_in_syndip/high_confidence_alleles_in_syndip:.1%}</td>""")
    table_html.append("</tr>")

    # row 2
    table_html.append("<tr>")
    table_html.append("""<td><i>2</i></td>""")
    table_html.append("""<td>Keep only insertions and deletions that represent <br />
    TR expansions or contractions based on these criteria:
    <ul>
        <li style="margin: 2px">span ≥ 9 base pairs</li>
        <li style="margin: 2px">have ≥ 3 repeats</li>
        <li style="margin: 2px">motif size ≥ 2bp and ≤ 50bp</li>
    </ul>
    </td>
    """)
    print(f"{format_n(len(df_variants))}    TOTAL STR loci")
    print(f"{format_n(len(df_alleles)) }    TOTAL STR alleles")

    table_html.append(f"""<td nowrap style="text-align: right">{len(df_variants):,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{len(df_alleles):,d}</td>""")
    table_html.append(f"""<td nowrap style="text-align: right">{len(df_alleles)/high_confidence_alleles_in_syndip:.1%}</td>""")
    table_html.append("</tr>")

    table_html.append("</table>")

    with open("figures_and_tables/table_template.html", "rt") as f, open(args.output_html, "wt") as fo:
        table_html_template = f.read()
        table_html_template = table_html_template % {
            "title": "Summary Stats",
            "table": "\n".join(table_html),
        }
        fo.write(table_html_template)


if __name__ == "__main__":
    main()
