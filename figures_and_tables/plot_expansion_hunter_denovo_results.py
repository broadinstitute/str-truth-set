import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "legend.fontsize": 13,
})

GREEN_COLOR = "#50AA44"


def bin_repeat_size_bp(long_allele_size_bp):
    limit = 500
    if long_allele_size_bp < 50:
        long_allele_size_bp = 50

    if long_allele_size_bp < limit:
        return str(25 * int(long_allele_size_bp/25))
    else:
        return f"{limit}+"


def compute_motif_size_bin(motif_size):
    if motif_size > 24:
        motif_size_bin = "25-50bp"
    elif motif_size > 6:
        motif_size_bin = "7-24bp"
    else:
        motif_size_bin = "2-6bp"

    return motif_size_bin


def compute_hue_column(row):
    if row["EHdn Concordance"] == "No Call":
        return "No matching EHdn call"

    if row["MotifSize"] > 24:
        motif_size_bin = "25-50bp"
    elif row["MotifSize"] > 6:
        motif_size_bin = "7-24bp"
    else:
        motif_size_bin = "2-6bp"

    return f"Matching EHdn call ({motif_size_bin} motif)"


def compute_concordance_summary(row):
    if row["EHdn Concordance With Truth Set"] == "Matching Expansion In Truth Set":
        return "Matches expansion in truth set"
    if not pd.isna(row["MatchingReferenceLocus"]):
        return "No expansion in truth set, but matches pure repeats in hg38"

    return "No expansion in truth set & no matching pure repeats in hg38"


def plot_truth_set_overlap_with_ehdn_calls(show_title=True):
    x_column = "RepeatSizeLongAllele (bin)"
    hue_column = "EHdn Concordance"

    truth_set_df = pd.read_table("~/code/str-truth-set/tool_comparison/results/expansion_hunter_denovo/"
                                 "CHM1_CHM13_WGS2.truth_set_EHdn_comparison_table.tsv")
    truth_set_df = truth_set_df[truth_set_df["IsPureRepeat"]]

    truth_set_df.loc[:, "MotifSize bin"] = truth_set_df.MotifSize.apply(compute_motif_size_bin)
    truth_set_df.loc[:, x_column] = truth_set_df["RepeatSizeLongAllele (bp)"].apply(bin_repeat_size_bp)
    truth_set_df.loc[:, hue_column] = truth_set_df.apply(compute_hue_column, axis=1)
    truth_set_df = truth_set_df.sort_values("RepeatSizeLongAllele (bp)", key=lambda values: [int(v) for v in values])

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8), dpi=80, tight_layout=True)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel("Truth Set Allele Size (bp)", labelpad=15, fontsize=13)
    ax.set_ylabel("Fraction of Truth Set Loci\nWith Matching ExpansionHunterDenovo Calls", labelpad=15, fontsize=13)

    sns.histplot(
        truth_set_df,
        x=x_column,
        hue=hue_column,
        binwidth=1,
        palette=["#c9c9c9", "#32CD32", "#299617", "#00703C"],
        multiple="fill", stat="proportion",
        discrete=True,
        legend=True,
        ax=ax)

    fig.tight_layout()

    if show_title:
        ax.set_title("Truth Set vs ExpansionHunterDenovo (EHdn) Calls", pad=75, fontsize=15)
    ax.get_legend().set_title(f"")

    # add # of loci label above each bar
    n_lookup = dict(truth_set_df.groupby(x_column).count().LocusId)
    for j, (xtick, text) in enumerate(zip(ax.get_xticks(), ax.get_xticklabels())):
        ax.text(xtick, 1.018, f"{n_lookup[text.get_text()]:,d} loci",
                ha="left", va="bottom", color="#777777", rotation=45, fontsize=12)

    ax.set_xticklabels([v if i%2 == 0 else "" for i, v in enumerate(ax.get_xticklabels())], fontsize=12)
    plt.yticks(fontsize=12)

    fig.tight_layout()

    ax.get_legend().set_bbox_to_anchor((0.51, 0.25))

    output_image_name = "truth_set_overlap_with_expansion_hunter_denovo_calls.svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_name}")


def plot_ehdn_calls_overlap_with_truth_set(show_title=True):
    ehdn_df = pd.read_table("~/code/str-truth-set/tool_comparison/results/expansion_hunter_denovo/"
                            "CHM1_CHM13_WGS2.EHdn_results_table.with_truth_set_concordance.tsv")

    ehdn_df.loc[:, "EHdn Concordance Summary"] = ehdn_df.apply(compute_concordance_summary, axis=1)

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 8), dpi=80, tight_layout=True)

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

    if show_title:
        ax.set_title("\n\nExpansionHunterDenovo Calls vs Truth Set", pad=20, fontsize=15)

    ax.get_legend().set_title(f"")
    ax.get_legend().set_frame_on(False)

    output_image_name = "expansion_hunter_denovo_calls_overlap_with_truth_set.svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_name}")


def main():
    plot_truth_set_overlap_with_ehdn_calls()
    plot_ehdn_calls_overlap_with_truth_set()


if __name__ == "__main__":
    main()