import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})


def plot_allele_size_distribution(df, plot_type=1, color_by=None, hue_order=None, is_pure_repeats=True):
    if is_pure_repeats:
        df = df[df.IsPureRepeat == "Yes"]
    else:
        df = df[df.IsPureRepeat == "No"]

    binwidth = None
    discrete = None
    bins = None
    if plot_type == 1:
        x_column = "NumRepeatsAwayFromReference"
        xlabel = "# of Repeats Relative to hg38"
        xlimit = 15
        minus_xlimit = -xlimit
        xticks = range(-xlimit, xlimit + 1, 2)
        xtick_labels = [f"{x}" if x < 0 else f"+{x}" for x in xticks]
        binwidth = 1
        discrete = True
        output_image_name = "allele_size_distribution_by_number_of_repeats"
        figure_title = "STR Allele Size Distribution"
        if color_by == "Multiallelic":
            figure_title = "Multiallelic Loci by Allele Size"
        elif color_by == "OverlapsSegDupIntervals":
            figure_title = "STR Overlap with Segmental Duplications"
    elif plot_type == 2:
        x_column = "NumBasePairsAwayFromReference"
        xlabel = "Allele Size (bp) Relative to hg38"
        xlimit = 39
        minus_xlimit = -xlimit
        xticks = range(-xlimit, xlimit + 1, 6)
        xtick_labels = [f"{x}" if x < 0 else f"+{x}" for x in xticks]
        bins = [b-1.5 for b in range(-xlimit, xlimit + 3, 3)]
        figure_title = "STR Allele Size Distribution"
        output_image_name = "allele_size_distribution_in_base_pairs"
    elif plot_type == 3:
        x_column = "MotifSize"
        xlabel = "Motif Size (bp)"
        xlimit = 50 if is_pure_repeats else 24
        minus_xlimit = 1
        xticks = [2, 3, 4, 5, 6] + list(range(9, xlimit + 1, 3))
        xtick_labels = xticks
        binwidth = 1
        discrete = True
        output_image_name = "allele_size_distribution_by_motif_size"
        figure_title = "STR Allele Size Distribution"
        if color_by == "Multiallelic":
            figure_title = "Multiallelic Loci by Motif Size"
    else:
        raise ValueError(f"Unexpected plot_type: {plot_type}")

    if is_pure_repeats:
        output_image_name += ".pure_repeats"
    else:
        output_image_name += ".with_interruptions"
        figure_title += "\n(interrupted repeats only)"

    _, ax = plt.subplots(figsize=(8, 7), dpi=80)
    ax.xaxis.labelpad = ax.yaxis.labelpad = 15
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Fraction of Alleles", fontsize=14)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize=13)
    ax.set_xlim(minus_xlimit - 0.52, xlimit + 0.52)

    plt.yticks(fontsize=13)
    p = sns.histplot(
        df,
        x=x_column,
        hue=color_by,
        hue_order=hue_order,
        bins=bins,
        binwidth=binwidth,
        multiple="stack" if not color_by else "fill",
        stat="proportion",
        discrete=discrete,
        ax=ax)

    p.set_title(figure_title, fontsize=14)

    if color_by:
        output_image_name += f".color_by_{color_by.lower()}"
        if color_by == "Multiallelic":
            ax.get_legend().set_title("Multiallelic")
        elif color_by == "OverlapsSegDupIntervals":
            ax.get_legend().set_title("\n".join([
                "STR Locus Overlaps ",
                "Segmental Duplication",
            ]))
        ax.get_legend().get_title().set_horizontalalignment('center')

        if plot_type == 3:
            if is_pure_repeats:
                ax.get_legend().set_bbox_to_anchor((0.9, 0.9))
            else:
                ax.get_legend().set_bbox_to_anchor((0.803, 0.9))

        output_image_name += ".svg"
        plt.savefig(f"{output_image_name}", bbox_inches="tight")
        plt.close()
        print(f"Saved {output_image_name}")

    print(f"Plotted {len(df):,d} allele records")


def plot_allele_size_distribution_x3(df, is_pure_repeats=True):
    """Plot allele size distribution split into 3 plots by motif size ranges"""
    if is_pure_repeats:
        df = df[df.IsPureRepeat == "Yes"]
    else:
        df = df[df.IsPureRepeat == "No"]

    fig, axes = plt.subplots(ncols=3, figsize=(30, 10), dpi=120, sharey=True)

    title = "Allele Size Distributions"
    if not is_pure_repeats:
        title += "\n(interrupted repeats only)"
    fig.suptitle(title, fontsize=24)

    for i, ax in enumerate(axes):

        xlimit = 1000
        step_size = 50

        if i == 0:
            df_current = df[(df.MotifSize >= 2) & (df.MotifSize <= 6)]
            label = "2bp to 6bp motifs"
        elif i == 1:
            df_current = df[(df.MotifSize >= 7) & (df.MotifSize <= 24)]
            label = "7bp to 24bp motifs"
        elif i == 2:
            label = "25bp to 50bp motifs"
            df_current = df[(df.MotifSize >= 25) & (df.MotifSize <= 50)]
        else:
            raise ValueError(f"{i}")

        ax.set_xlabel("Allele Size (bp)", labelpad=15, fontsize=18)
        if i == 0:
            ax.set_ylabel("Number of Alleles", labelpad=15, fontsize=18)
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xticks(range(-xlimit, xlimit + 1, step_size*4))
        ax.set_xticklabels([f"{x}" if x < 0 else f"+{x}" for x in range(-xlimit, xlimit + 1, step_size*4)])
        ax.set_xlim(-xlimit - 0.52, xlimit + 0.52)
        ax.set_xticklabels(
            ax.get_xticklabels(),
            rotation=45,
            fontsize=18,
            horizontalalignment="right",
            rotation_mode='anchor')
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=18)
        ax.yaxis.set_tick_params(labelleft=True)

        p = sns.histplot(
            df_current,
            x="NumBasePairsAwayFromReference",
            bins=[b-1.5 for b in range(-xlimit, xlimit + step_size, step_size)],
            ax=ax)

        p.set(yscale='log')
        p.set_title(f"{len(df_current):,d} alleles at {len(set(df_current.LocusId)):,d} loci ({label})", fontsize=18)

        output_image_name = "allele_size_distribution_by_number_of_repeats.x3"
        if is_pure_repeats:
            output_image_name += ".pure_repeats"
        else:
            output_image_name += ".with_interruptions"
        output_image_name += ".svg"

        fig.savefig(f"{output_image_name}", bbox_inches="tight")
        print(f"Saved {output_image_name}")

        print(f"Plotted {len(df):,d} allele records")


def plot_motif_distribution(df, is_pure_repeats=True):
    print("Plotting allele distribution by motif size")
    if is_pure_repeats:
        df = df[df.IsPureRepeat == "Yes"]
    else:
        df = df[df.IsPureRepeat == "No"]

    hue_limit = 5
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReference"]
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = np.where(
        df["NumRepeatsAwayFromReferenceTruncated"] > hue_limit, hue_limit, df["NumRepeatsAwayFromReferenceTruncated"])
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = np.where(
        df["NumRepeatsAwayFromReferenceTruncated"] < -hue_limit, -hue_limit, df["NumRepeatsAwayFromReferenceTruncated"])
    df = df.sort_values("NumRepeatsAwayFromReferenceTruncated")

    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].astype("int").apply(
        lambda i: f"{i}" if i < 0 else f"+{i}")
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].replace(
        f"-{hue_limit}", f"-{hue_limit} or less")
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].replace(
        f"+{hue_limit}", f"+{hue_limit} or more")

    _, axes = plt.subplots(2, 1, figsize=(8, 12), dpi=80)

    title = "# of Repeats in hg38 at Truth Set STR Loci"
    if not is_pure_repeats:
        title += "\n(interrupted repeats only)"
    ax = axes[0]
    p = sns.histplot(
        df,
        x="NumRepeatsInReference",
        binwidth=1,
        discrete=True,
        ax=ax)

    p.set_title(title, fontsize=14)
    ax.set_xlabel("", fontsize=14)
    ax.set_yscale("log")
    ax.set_ylabel("# of Loci", fontsize=14)

    ax = axes[1]
    sns.histplot(
        df,
        x="NumRepeatsInReference",
        hue="NumRepeatsAwayFromReferenceTruncated",
        palette=sns.color_palette("bwr_r", 10)[::-1],
        binwidth=1,
        discrete=True,
        stat="proportion",
        multiple="fill",
        ax=ax)

    xlimit = 51
    for ax in axes:
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xticks(range(0, xlimit + 1, 5))
        ax.set_xlim(-1, xlimit)

    ax.set_xlabel("# of Repeats in hg38", fontsize=14)
    ax.set_ylabel("Fraction of Alleles", fontsize=14)
    plt.rcParams.update({"text.usetex": True})
    ax.get_legend().set_title("\n".join([
        "# of Repeats in",
        "Truth Set Allele",
        "Relative to hg38",
        "$_{contractions\ are < 0}$",
        "$_{expansions\ are > 0}   $",
    ]), prop={'size': 12})
    ax.get_legend()._legend_box.align = "left"
    ax.get_legend().set_bbox_to_anchor((1.35, 0.25))
    ax.get_legend().get_frame().set_color("white")

    output_image_name = "reference_locus_size_distribution"
    if is_pure_repeats:
        output_image_name += ".pure_repeats"
    else:
        output_image_name += ".with_interruptions"
    output_image_name += ".svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_name}")

    print(f"Plotted {len(df):,d} allele records")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--skip-plot1", action="store_true")
    p.add_argument("--skip-plot2", action="store_true")
    args = p.parse_args()

    input_table_path = "../STR_truthset.v1.alleles.tsv.gz"

    print(f"Reading {input_table_path}")
    df = pd.read_table(input_table_path)
    print(f"Parsed {len(df):,d} rows from {input_table_path}")
    df = df.rename(columns={"IsMultiallelic": "Multiallelic"})
    df.loc[:, "NumRepeatsAwayFromReference"] = df.NumRepeats - df.NumRepeatsInReference
    df.loc[:, "NumBasePairsAwayFromReference"] = df.NumRepeatsAwayFromReference * df.MotifSize

    if not args.skip_plot1:
        plot_allele_size_distribution(df, plot_type=1, is_pure_repeats=True)
        plot_allele_size_distribution(df, plot_type=1, is_pure_repeats=False)
        plot_allele_size_distribution(df, plot_type=2, is_pure_repeats=True)
        plot_allele_size_distribution(df, plot_type=2, is_pure_repeats=False)
        plot_allele_size_distribution(df, plot_type=1, color_by="Multiallelic", hue_order=["No", "Yes"], is_pure_repeats=True)
        plot_allele_size_distribution(df, plot_type=1, color_by="Multiallelic", hue_order=["No", "Yes"], is_pure_repeats=False)
        plot_allele_size_distribution(df, plot_type=3, color_by="Multiallelic", hue_order=["No", "Yes"], is_pure_repeats=True)
        plot_allele_size_distribution(df, plot_type=3, color_by="Multiallelic", hue_order=["No", "Yes"], is_pure_repeats=False)

        plot_allele_size_distribution(df, plot_type=1, color_by="OverlapsSegDupIntervals", hue_order=["No", "Yes"], is_pure_repeats=True)
        plot_allele_size_distribution_x3(df, is_pure_repeats=True)

    if not args.skip_plot2:
        plot_motif_distribution(df, is_pure_repeats=True)


if __name__ == "__main__":
    main()