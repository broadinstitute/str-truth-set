import argparse
from numbers_utils import search, format_n, format_np
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants-table", help="Path of variants table",  default="STR_truth_set.v1.variants.tsv.gz")
    parser.add_argument("--step-A-log",     help="Path of step A log file", default="step_A.log")
    parser.add_argument("--output-html",   help="Path of output table",    default="figures_and_tables/table3_repeat_catalogs.html")
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


    print("-"*100)
    print("Numbers for intro:")
    print("-"*100)
    total = len(df_variants[df_variants.IsFoundInReference])
    gangstr_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from GangSTRCatalog17", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    illumina_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from IlluminaSTRCatalog", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    hipstr_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from HipSTRCatalog", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_6bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats6bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_9bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats9bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_12bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats12bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)

    missed_by_gangstr_catalog = search(
        f"step8:STR:[ ]+GangSTRCatalog17.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_gangstr_catalog, total)} missed by GangSTRCatalog17 catalog")

    missed_by_illumina_catalog = search(
        f"step8:STR:[ ]+IlluminaSTRCatalog.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_illumina_catalog, total)} missed by IlluminaSTRCatalog catalog")

    missed_by_hipstr_catalog = search(
        f"step8:STR:[ ]+HipSTRCatalog.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_hipstr_catalog, total)} missed by HipSTRCatalog catalog")

    missed_by_6bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats6bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_6bp_catalog, total)} missed by TRFPureRepeats6bp catalog")

    missed_by_9bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats9bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_9bp_catalog, total)} missed by TRFPureRepeats9bp catalog")

    missed_by_12bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats12bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_12bp_catalog, total)} missed by TRFPureRepeats12bp catalog")

    header = [
        #"",
        "Catalog name",
        "How it was created",
        "Catalog size<br/>(# of loci)",
        "# of truth set<br/>loci missed",
        "% of truth set<br />loci missed",
    ]

    data = [
        [
            "GangSTR catalog v17",
            "Running TandemRepeatFinder (TRF) on hg38 <br />"
            "with mismatch penalty = 5, indel penalty = 17<br />",
            gangstr_catalog_size,
            missed_by_gangstr_catalog,
            missed_by_gangstr_catalog/total,
        ],
        [
            "Illumina catalog",
            "polymorphic STR loci based on <br />"
            "2,504 genomes from 1kGP",
            illumina_catalog_size,
            missed_by_illumina_catalog,
            missed_by_illumina_catalog/total,
        ],
        [
            "HipSTR catalog",
            "Running TRF on hg38 using default params:<br />"
            "mismatch penalty = 7, indel penalty = 7",
            hipstr_catalog_size,
            missed_by_hipstr_catalog,
            missed_by_hipstr_catalog/total,
        ],
        [
            "New catalog (large)",
            "Running TRF on hg38 using very large<br />"
            "indel & mismatch penalties to find<br />"
            "all pure repeats that <b>span at least 6bp</b>",
            new_6bp_catalog_size,
            missed_by_6bp_catalog,
            missed_by_6bp_catalog/total,
        ],
        [
            "New catalog",
            "Running TRF on hg38 using very large<br />"
            "indel & mismatch penalties to find<br />"
            "all pure repeats that <b>span at least 9bp</b>",
            new_9bp_catalog_size,
            missed_by_9bp_catalog,
            missed_by_9bp_catalog/total,
        ],
        [
            "New catalog (smaller)",
            "Running TRF on hg38 using very large<br />"
            "indel & mismatch penalties to find<br />"
            "all pure repeats that <b>span at least 12bp</b>",
            new_12bp_catalog_size,
            missed_by_12bp_catalog,
            missed_by_12bp_catalog/total,
        ],
    ]

    table_html = []
    table_html.append(f"<table>")
    table_html.append(f"""<tr>{''.join(f'<th nowrap style="text-align: center">{h}</th>' for h in header)}</tr>""")

    for i, data_row in enumerate(data):
        table_html.append(f"<tr>"
                          #f"""<td style="line-height: 1.5; vertical-align: top">{i+1}</td>"""
                          f"""<td style="line-height: 1.5; vertical-align: top">{data_row[0]}</td>"""
                          f"""<td style="line-height: 1.5; vertical-align: top">{data_row[1]}</td>"""
                          f"""<td style="text-align: right; line-height: 1.5; vertical-align: top">{data_row[2]:,d}</td>"""
                          f"""<td style="text-align: right; line-height: 1.5; vertical-align: top">{data_row[3]:,d}</td>"""
                          f"""<td style="text-align: right; line-height: 1.5; vertical-align: top">{data_row[4]:0.1%}</td>"""
                          f"</tr>")

    table_html.append("</table>")

    with open("figures_and_tables/table_template.html", "rt") as f, open(args.output_html, "wt") as fo:
        table_html_template = f.read()
        table_html_template = table_html_template % {
            "title": "Completeness of available TR catalogs",
            "table": "\n".join(table_html),
        }
        fo.write(table_html_template)


if __name__ == "__main__":
    main()