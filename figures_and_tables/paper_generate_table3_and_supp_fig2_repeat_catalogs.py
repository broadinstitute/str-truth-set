import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pybedtools
import seaborn as sns

from numbers_utils import search, format_np

new_catalogs_name = "repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_%dbp.bed.gz"


def compute_catalog_comparison_table(args):
    total_loci_in_truth_set = 0
    truth_set_tmp_path = os.path.join(args.temp_dir, "truth_set_loci.present_in_hg38.bed")
    with gzip.open(args.truth_set_variants_bed, "rt") as f, open(truth_set_tmp_path, "wt") as f_out:
        for line in f:
            total_loci_in_truth_set += 1
            fields = line.strip().split("\t")
            start_0based = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            if end - start_0based > 0:
                # this means the locus is not novel - ie. is present in the reference genome
                f_out.write(line)
            else:
                # for truth set loci that are novel, the summary string (stored in the name field of the bed file) should
                # look like "2bp:TG:INS:0=>11:HET". As a conistency check, make sure it includes the "0=>" string indicating
                # it has 0 repeats in the reference genome.
                if "0=>" not in name:
                    print("ERROR: unexpected summary string for non-ref truth set locus:", fields)

    print(f"Read {total_loci_in_truth_set:,d} loci from {args.truth_set_variants_bed} that are present in hg38")

    catalog_paths = [
        ("Illumina catalog", args.illumina_catalog_bed, None),
        ("GangSTR v17 catalog", args.gangstr_v17_catalog_bed, None),
        ("HipSTR catalog", args.hipstr_catalog_bed, None),
    ]

    for locus_spans_min_base_pairs in range(6, 25, 3):
        new_catalog_path = os.path.join(args.new_catalogs_dir, new_catalogs_name % locus_spans_min_base_pairs)
        if not os.path.isfile(new_catalog_path):
            print("ERROR: new catalog file not found:", new_catalog_path,
                  "Run ./ref/other/generate_repeat_catalogs_spanning_x_bp.sh to generate it.")
            continue
        catalog_paths.append((
            f"New catalog",
            new_catalog_path,
            locus_spans_min_base_pairs,
        ))

    table_rows = []
    for label, catalog_path, locus_spans_min_base_pairs in catalog_paths:
        catalog_size = pybedtools.BedTool(catalog_path).count()

        # use pybedtools to see how many truth set loci don't overlap any loci in the catalog.
        #   The default minimum overlap is 1bp.
        #   A=True means "Remove entire feature if any overlap."
        num_loci_missed_by_catalog = pybedtools.BedTool(truth_set_tmp_path).subtract(catalog_path, A=True).count()

        print(f"{num_loci_missed_by_catalog:,d} out of {total_loci_in_truth_set:,d} ({num_loci_missed_by_catalog / total_loci_in_truth_set:.1%}) "
              f"loci missed by {label} which contains {catalog_size:,d} loci")

        table_rows.append({
            "catalog_long_name": label + (" (loci spanning at least %dbp)" % locus_spans_min_base_pairs if locus_spans_min_base_pairs else ""),
            "catalog": label,
            "catalog_size": catalog_size,
            "num_loci_missed_by_catalog": num_loci_missed_by_catalog,
            "fraction_missed_by_catalog": num_loci_missed_by_catalog / total_loci_in_truth_set,
            "truth_set_size": total_loci_in_truth_set,
            "locus_spans_min_base_pairs": locus_spans_min_base_pairs,
        })

    df = pd.DataFrame(table_rows)
    return df


def generate_catalog_comparison_plot(df, args):
    df = df[df['locus_spans_min_base_pairs'].isna() | (df['locus_spans_min_base_pairs'] <= 24)]

    df["catalog"] = df["catalog"].replace("New catalog", "New catalogs")

    # use seaborn to draw a scatter plot of the catalog size vs. the fraction of truth set loci missed by the catalog
    # the points should be labeled with the catalog name
    # the x-axis should be labeled "catalog size" and should have a log scale
    # the y-axis should be labeled "fraction of truth set loci missed by catalog"
    # the plot should be saved to "catalog_comparison.png"

    sns.set_style("whitegrid")
    sns.set_context("paper", font_scale=1.5)
    fig, ax = plt.subplots(figsize=(10, 6), tight_layout=True)
    sns.scatterplot(
        data=df,
        x="catalog_size",
        y="fraction_missed_by_catalog",
        hue="catalog",
        # set the dot size to larger
        s=100,
        linewidth=2,
        ax=ax)

    ax.set_xscale("log")
    ax.set_xticks([1e5, 1e6, 1e7])
    ax.set_xlim(7e4, 1e7)
    y_limit = 0.5
    ax.set_ylim(0, y_limit)

    ax.set_yticks(np.arange(0, y_limit + 0.1, 0.1))
    ax.tick_params(axis='both', which='major', labelsize=16)

    # label the illumina catalog point with the catalog name
    for _, row in df.iterrows():
        is_new_catalog = row["catalog"].lower().startswith("new")
        label_y = row["fraction_missed_by_catalog"]
        if is_new_catalog:
            label = f"≥ {int(row['locus_spans_min_base_pairs'])}bp "
            #if row["locus_spans_min_base_pairs"] == 24:
            #    label = "Loci " + label
            label_y -= 0.025
        else:
            label = f" {row['catalog']}"
            label_y += 0.015

        ax.text(
            row["catalog_size"],
            label_y,
            label,
            horizontalalignment='right' if is_new_catalog else 'left',
            fontsize=14)
        #ax.text(row["catalog_size"], row["fraction_missed_by_catalog"], row["catalog"], fontsize=12)

    ax.xaxis.labelpad = 15
    ax.yaxis.labelpad = 15

    ax.legend(loc='center left', handlelength=1.2, bbox_to_anchor=(1, 0.88), fontsize=14)
    ax.get_legend().set_title(None)
    ax.get_legend().set_frame_on(False)

    ax.set_xlabel("Catalog Size (# of loci)", fontsize=16)
    ax.set_ylabel("Fraction of Truth Set Loci Missed by Catalog", fontsize=16)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(color='grey', linestyle='-', linewidth=0.25, alpha=0.5)

    output_path = os.path.join(args.output_dir, f"catalog_comparison.{args.image_type}")
    plt.savefig(output_path)
    plt.close()

    print(f"Saved plot to {output_path}")


def generate_html_table(df, args):
    df = df.set_index("catalog_long_name")

    header = [
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
            df.loc["Illumina catalog", "catalog_size"],
            df.loc["Illumina catalog", "num_loci_missed_by_catalog"],
            df.loc["Illumina catalog", "fraction_missed_by_catalog"],
        ],
        [
            "Illumina catalog",
            "polymorphic STR loci based on <br />"
            "2,504 genomes from 1kGP",
            df.loc["GangSTR v17 catalog", "catalog_size"],
            df.loc["GangSTR v17 catalog", "num_loci_missed_by_catalog"],
            df.loc["GangSTR v17 catalog", "fraction_missed_by_catalog"],
        ],
        [
            "HipSTR catalog",
            "Running TRF on hg38 using default params:<br />"
            "mismatch penalty = 7, indel penalty = 7",
            df.loc["HipSTR catalog", "catalog_size"],
            df.loc["HipSTR catalog", "num_loci_missed_by_catalog"],
            df.loc["HipSTR catalog", "fraction_missed_by_catalog"],
        ],
    ]
    for min_bp in [6, 9, 12, 15]:
        data.append([
            f"New catalog<br />(loci ≥ {min_bp}bp)",
            "Running TRF on hg38 using very large<br />"
            "indel & mismatch penalties to find<br />"
            f"all pure repeats that <b>span at least {min_bp}bp</b>",
            df.loc[f"New catalog (loci spanning at least {min_bp}bp)", "catalog_size"],
            df.loc[f"New catalog (loci spanning at least {min_bp}bp)", "num_loci_missed_by_catalog"],
            df.loc[f"New catalog (loci spanning at least {min_bp}bp)", "fraction_missed_by_catalog"],
        ])


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

    with open(args.table_template_html, "rt") as f, open(os.path.join(args.output_dir, args.output_html), "wt") as fo:
        table_html_template = f.read()
        table_html_template = table_html_template % {
            "title": "Completeness of available TR catalogs",
            "table": "\n".join(table_html),
        }
        fo.write(table_html_template)

    print("Wrote html table to", os.path.join(args.output_dir, args.output_html))


def check_numbers(df, args):
    """As a sanity check, compare numbers to those in step_A.log"""

    with open(args.step_A_log, "rt") as f:
        stepA_log_contents = f.read()

    df_variants = pd.read_table(args.variants_table)

    total = len(df_variants[df_variants.IsFoundInReference])
    gangstr_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from GangSTRCatalog17", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    illumina_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from IlluminaSTRCatalog", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    hipstr_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from HipSTRCatalog", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_6bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats6bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_9bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats9bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)
    new_12bp_catalog_size = search("step8:STR:[ ]+Parsed[ ]+([0-9,]+) intervals from TRFPureRepeats12bp", stepA_log_contents, type=int, expected_number_of_matches=2, use_match_i=0)

    df = df.set_index("catalog_long_name")


    print("-"*100)
    print("Checking numbers vs those in step_A.log:")
    print("-"*100)

    missed_by_gangstr_catalog = search(
        f"step8:STR:[ ]+GangSTRCatalog17.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_gangstr_catalog, total)} missed by GangSTRCatalog17 catalog")
    if df.loc["GangSTR v17 catalog", "catalog_size"] != gangstr_catalog_size:
        raise ValueError(f"GangSTR catalog size mismatch: {df.loc['GangSTR v17 catalog', 'catalog_size']} != {gangstr_catalog_size}")

    missed_by_illumina_catalog = search(
        f"step8:STR:[ ]+IlluminaSTRCatalog.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_illumina_catalog, total)} missed by IlluminaSTRCatalog catalog")
    if df.loc["Illumina catalog", "catalog_size"] != illumina_catalog_size:
        raise ValueError(f"Illumina catalog size mismatch: {df.loc['Illumina catalog', 'catalog_size']} != {illumina_catalog_size}")

    missed_by_hipstr_catalog = search(
        f"step8:STR:[ ]+HipSTRCatalog.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_hipstr_catalog, total)} missed by HipSTRCatalog catalog")
    if df.loc["HipSTR catalog", "catalog_size"] != hipstr_catalog_size:
        raise ValueError(f"HipSTR catalog size mismatch: {df.loc['HipSTR catalog', 'catalog_size']} != {hipstr_catalog_size}")

    missed_by_6bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats6bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_6bp_catalog, total)} missed by TRFPureRepeats6bp catalog")
    key = "New catalog (loci spanning at least 6bp)"
    if df.loc[key, "catalog_size"] != new_6bp_catalog_size:
        raise ValueError(f"TRF 6bp catalog size mismatch: {df.loc[key, 'catalog_size']} != {new_6bp_catalog_size}")

    missed_by_9bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats9bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_9bp_catalog, total)} missed by TRFPureRepeats9bp catalog")
    key = f"New catalog (loci spanning at least 9bp)"
    if df.loc[key, "catalog_size"] != new_9bp_catalog_size:
        raise ValueError(f"TRF 9bp catalog size mismatch: {df.loc[key, 'catalog_size']} != {new_9bp_catalog_size}")

    missed_by_12bp_catalog = search(
        f"step8:STR:[ ]+TRFPureRepeats12bp.*[ ]([0-9,]+) out of [ ]+{total:,d} [(][ ]+[0-9]+[.][0-9]+[%][)].* truth set loci:.*no overlap",
        stepA_log_contents, type=int)
    print(f"{format_np(missed_by_12bp_catalog, total)} missed by TRFPureRepeats12bp catalog")
    key = f"New catalog (loci spanning at least 12bp)"
    if df.loc[key, "catalog_size"] != new_12bp_catalog_size:
        raise ValueError(f"TRF 12bp catalog size mismatch: {df.loc[key, 'catalog_size']} != {new_12bp_catalog_size}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step-A-log", default="../step_A.log")
    parser.add_argument("--variants-table", help="Path of variants table",  default="../STR_truth_set.v1.variants.tsv.gz")
    parser.add_argument("--truth-set-variants-bed", default="../STR_truth_set.v1.variants.bed.gz")
    parser.add_argument("--illumina-catalog-bed", default="../ref/other/illumina_variant_catalog.sorted.bed.gz")
    parser.add_argument("--gangstr-v17-catalog-bed", default="../ref/other/hg38_ver17.adjusted.bed.gz")
    parser.add_argument("--hipstr-catalog-bed", default="../ref/other/hg38.hipstr_reference.adjusted.bed.gz")
    parser.add_argument("--new-catalogs-dir", default="../ref/other/")
    parser.add_argument("--table-template-html", default="table_template.html")
    parser.add_argument("--temp-dir", default="/tmp")
    parser.add_argument("--image-type", default="png", choices=["svg", "png"])
    parser.add_argument("--output-dir", default=".")
    parser.add_argument("--output-html",   help="Path of output table",    default="table3_repeat_catalogs.html")

    args = parser.parse_args()
    df = compute_catalog_comparison_table(args)
    check_numbers(df, args)
    generate_html_table(df, args)
    generate_catalog_comparison_plot(df, args)


if __name__ == "__main__":
    main()