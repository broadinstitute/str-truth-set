import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
from pprint import pprint
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
    "legend.fontsize": 13,
})


def bin_repeat_size_bp(long_allele_size_bp):
    limit = 500
    if long_allele_size_bp < 50:
        long_allele_size_bp = 50

    if long_allele_size_bp < limit:
        return str(25 * int(long_allele_size_bp/25))
    else:
        return f"{limit}+"


def compute_motif_size_bin(motif_size):
    if motif_size >= 25:
        motif_size_bin = "25-50bp"
    elif motif_size >= 7:
        motif_size_bin = "7-24bp"
    elif motif_size >= 2:
        motif_size_bin = "2-6bp"
    else:
        motif_size_bin = "1bp"

    return motif_size_bin


def compute_hue_column(row):
    if row["EHdn Concordance"] == "No Call":
        prefix = "No matching EHdn call"
    else:
        prefix = "Matching EHdn call"

    if row["MotifSize"] >= 25:
        motif_size_bin = "25-50bp"
    elif row["MotifSize"] >= 7:
        motif_size_bin = "7-24bp"
    elif row["MotifSize"] >= 2:
        motif_size_bin = "2-6bp"
    else:
        motif_size_bin = "1bp"

    #if row["EHdn Concordance"] == "No Call" and row["RepeatSizeLongAllele (bp)"] >= 250 and row["MotifSize"] < 25:
    #    pprint(row.to_dict())
    if row["EHdn Concordance"] != "No Call" and row["RepeatSizeLongAllele (bp)"] >= 250 and row["MotifSize"] <= 6:
        pprint(row.to_dict())

    return f"{prefix} ({motif_size_bin} motif)"


def compute_concordance_summary(row):
    if row["EHdn Concordance With Truth Set"] == "Matching Expansion In Truth Set":
        return "Matches expansion in truth set"
    if not pd.isna(row["MatchingReferenceLocus"]):
        return "No expansion in truth set, but matches pure repeats in hg38"

    return "No expansion in truth set & no matching pure repeats in hg38"


def plot_truth_set_overlap_with_ehdn_calls(args):
    x_column = "RepeatSizeLongAllele (bin)"
    hue_column = "EHdn Concordance (bin)"

    filename_suffix = ""

    truth_set_df = pd.read_table(args.truth_set_ehdn_input_table)
    if args.only_pure_repeats:
        filename_suffix += ".pure_repeats"
        truth_set_df = truth_set_df[truth_set_df["IsPureRepeat"]]

    truth_set_df.loc[:, "MotifSize bin"] = truth_set_df.MotifSize.apply(compute_motif_size_bin)
    truth_set_df.loc[:, x_column] = truth_set_df["RepeatSizeLongAllele (bp)"].apply(bin_repeat_size_bp)
    truth_set_df.loc[:, hue_column] = truth_set_df.apply(compute_hue_column, axis=1)
    truth_set_df = truth_set_df.sort_values("RepeatSizeLongAllele (bp)")

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(args.width, args.height), dpi=80, tight_layout=True)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel("Truth Set Total Allele Size (bp)", labelpad=15, fontsize=13)
    ax.set_ylabel("Fraction of Truth Set Expansions\nWith Matching ExpansionHunterDenovo Calls", labelpad=15, fontsize=13)

    sns.histplot(
        truth_set_df,
        x=x_column,
        hue=hue_column,
        hue_order=[
            "No matching EHdn call (2-6bp motif)",
            "No matching EHdn call (7-24bp motif)",
            "No matching EHdn call (25-50bp motif)",
            "Matching EHdn call (2-6bp motif)",
            "Matching EHdn call (7-24bp motif)",
            "Matching EHdn call (25-50bp motif)",
        ],
        binwidth=1,
        palette=[
            "#4C698E",
            "#72A1D1",
            "#c9c9c9",
            "#32CD32",
            "#299617",
            "#00703C",
        ],
        multiple="fill", stat="proportion",
        discrete=True,
        legend=True,
        ax=ax)

    fig.tight_layout()

    if args.show_title:
        ax.set_title("Truth Set vs ExpansionHunterDenovo (EHdn) Calls", pad=75, fontsize=15)
    ax.get_legend().set_title(f"")

    # add # of loci label above each bar
    n_lookup = dict(truth_set_df.groupby(x_column).count().LocusId)
    for j, (xtick, text) in enumerate(zip(ax.get_xticks(), ax.get_xticklabels())):
        ax.text(xtick, 1.018, f"{n_lookup[text.get_text()]:,d}",
                ha="left", va="bottom", color="#777777", rotation=45, fontsize=12)

    ax.text(1, 1.1, "Truth Set Expansion Variants Per Bin",
            ha="right", va="bottom", color="#777777", fontsize=12, transform=ax.transAxes)

    ax.set_xticklabels([v if i%2 == 0 else "" for i, v in enumerate(ax.get_xticklabels())], fontsize=12)
    plt.yticks(fontsize=12)

    fig.tight_layout()

    sns.move_legend(ax, loc="lower right")

    output_path = os.path.join(args.output_dir, f"truth_set_overlap_with_expansion_hunter_denovo_calls{filename_suffix}.{args.image_type}")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Plotted {len(truth_set_df):,d} truth set alleles")
    print(f"Saved {output_path}")


def plot_ehdn_calls_overlap_with_truth_set(args):
    ehdn_df = pd.read_table(args.ehdn_truth_set_input_table)

    ehdn_df.loc[:, "EHdn Concordance Summary"] = ehdn_df.apply(compute_concordance_summary, axis=1)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(args.width, args.height), tight_layout=True)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel("Motif Size (bp)", labelpad=15, fontsize=14)
    ax.set_ylabel("ExpansionHunterDenovo Call Counts", labelpad=15, fontsize=14)
    xticks = [2, 3, 4, 5, 6] + list(range(10, 51, 5))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=12)
    plt.yticks(fontsize=12)
    sns.histplot(
        ehdn_df,
        x="MotifSize",
        hue="EHdn Concordance Summary",
        binwidth=1,
        palette=["#FF5347", "#FFC939", "#689F38"],
        multiple="stack",
        stat="count",
        discrete=True,
        legend=True,
        ax=ax)

    if args.show_title:
        ax.set_title("\n\nExpansionHunterDenovo Calls vs Truth Set", pad=20, fontsize=15)

    ax.get_legend().set_title(f"")
    ax.get_legend().set_frame_on(False)

    output_path = os.path.join(args.output_dir, f"expansion_hunter_denovo_calls_overlap_with_truth_set.{args.image_type}")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()

    print(f"Saved {output_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default=".")
    p.add_argument("--width", default=10, type=float)
    p.add_argument("--height", default=8, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"], help="Image type to generate")
    p.add_argument("--show-title", action="store_true", help="Show plot title")

    g = p.add_argument_group("Filters")
    g.add_argument("--only-pure-repeats", action="store_true", help="Plot only loci with pure repeats")

    p.add_argument("--truth-set-ehdn-input-table",
                   default="../tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/CHM1_CHM13_WGS2.downsampled_to_30x.truth_set_EHdn_comparison_table.tsv")
    p.add_argument("--ehdn-truth-set-input-table",
                   default="../tool_comparison/results_for_downsampled_30x_bam/expansion_hunter_denovo/CHM1_CHM13_WGS2.downsampled_to_30x.EHdn_results_table.with_truth_set_concordance.tsv")
    args = p.parse_args()

    plot_truth_set_overlap_with_ehdn_calls(args)
    plot_ehdn_calls_overlap_with_truth_set(args)


if __name__ == "__main__":
    main()