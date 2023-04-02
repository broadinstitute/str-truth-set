import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
})


def plot_allele_size_distribution(df_truth_set, args, plot_type=1, color_by=None, hue_order=None, yaxis_fraction=False, palette=None):
    binwidth = None
    discrete = None
    bins = None
    if plot_type == 1:
        x_column = "NumRepeatsAwayFromReference"
        xlabel = "Allele Size (# of Repeats)"
        xlimit = 15
        minus_xlimit = -xlimit
        xticks = range(-xlimit, xlimit + 1, 2)
        xtick_labels = [f"{x}" if x < 0 else f"+{x}" for x in xticks]
        binwidth = 1
        discrete = True
        output_image_name = "allele_size_distribution_by_number_of_repeats"
        figure_title = "Allele Size Distribution"
        if color_by == "Multiallelic":
            figure_title = "Multiallelic Loci by Allele Size"
        elif color_by == "OverlapsSegDupIntervals":
            figure_title = "Overlap with Segmental Duplications"
    elif plot_type == 2:
        x_column = "NumBasePairsAwayFromReference"
        xlabel = "Allele Size (bp)"
        xlimit = 39
        minus_xlimit = -xlimit
        xticks = range(-xlimit, xlimit + 1, 6)
        xtick_labels = [f"{x}" if x < 0 else f"+{x}" for x in xticks]
        bins = [b-1.5 for b in range(-xlimit, xlimit + 3, 3)]
        figure_title = "Allele Size Distribution"
        output_image_name = "allele_size_distribution_in_base_pairs"
    elif plot_type == 3:
        x_column = "MotifSize"
        xlabel = "Motif Size (bp)"
        xlimit = 24
        minus_xlimit = 1
        xticks = [2, 3, 4, 5, 6] + list(range(9, xlimit + 1, 3))
        xtick_labels = xticks
        binwidth = 1
        discrete = True
        output_image_name = "allele_size_distribution_by_motif_size"
        figure_title = "Allele Size Distribution"
        if color_by == "Multiallelic":
            figure_title = "Multiallelic Loci by Motif Size"
    else:
        raise ValueError(f"Unexpected plot_type: {plot_type}")

    if args.only_pure_repeats:
        output_image_name += ".only_pure_repeats"

    _, ax = plt.subplots(figsize=(8, 7) if args.width is None or args.height is None else (args.width, args.height))
    ax.xaxis.labelpad = ax.yaxis.labelpad = 15
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel("Fraction of Alleles", fontsize=14)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels, fontsize=12)
    ax.set_xlim(minus_xlimit - 0.52, xlimit + 0.52)

    plt.yticks(fontsize=13)
    p = sns.histplot(
        df_truth_set,
        x=x_column,
        hue=color_by,
        hue_order=hue_order,
        bins=bins,
        binwidth=binwidth,
        multiple="fill" if yaxis_fraction else "stack",
        stat="proportion",
        discrete=discrete,
        palette=palette,
        ax=ax)


    if args.show_title:
        p.set_title(figure_title, fontsize=14)

    if color_by:
        output_image_name += f".color_by_{color_by.lower().replace(' ', '_')}"
        if color_by == "Multiallelic":
            ax.get_legend().set_title("Multiallelic")
        elif color_by == "OverlapsSegDupIntervals":
            ax.get_legend().set_title("\n".join([
                "TR Locus Overlaps ",
                "Segmental Duplication",
            ]))
        ax.get_legend().get_title().set_horizontalalignment('center')
        if plot_type == 1:
            sns.move_legend(ax, loc="upper left")
        elif plot_type == 3:
            sns.move_legend(ax, loc="upper right")
        ax.get_legend().set_frame_on(False)

    output_path = os.path.join(args.output_dir, output_image_name + f".{args.image_type}")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {output_path}")

    print(f"Plotted {len(df_truth_set):,d} allele records")


def plot_allele_size_and_motif_distribution(df_truth_set, args, color_by=None, hue_order=None, palette=None):
    output_image_name = "allele_size_distribution_by_number_of_repeats_and_motif_size"
    if args.only_pure_repeats:
        output_image_name += ".only_pure_repeats"

    figure_title = "Allele Size Distribution"
    if color_by == "Multiallelic":
        figure_title = "Multiallelic Loci"
    elif color_by == "OverlapsSegDupIntervals":
        figure_title = "Overlap with Segmental Duplications"

    fig, axes = plt.subplots(1, 2, figsize=(16, 6) if args.width is None or args.height is None else (args.width, args.height))
    if args.show_title:
        suptitle_artist = fig.suptitle(figure_title, fontsize=15)
        extra_artists = [suptitle_artist]
    else:
        extra_artists = []

    for i, ax in enumerate(axes):
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        if i == 0:
            xlimit = 15
            ax.set_xlabel("Allele Size (# of Repeats)", fontsize=13)
            ax.set_xticks(range(-xlimit, xlimit + 1, 2))
            ax.set_xticklabels([f"{x}" if x < 0 else f"+{x}" for x in range(-xlimit, xlimit + 1, 2)], fontsize=12)
            ax.set_xlim(-xlimit - 0.52, xlimit + 0.52)
        else:
            xlimit = 24
            ax.set_xlabel("Motif Size (bp)", fontsize=13)
            ax.set_xticks(range(2, xlimit + 1, 1))
            ax.set_xticklabels([f"{x}" if x <= 6 or x % 3 == 0 else "" for x in range(2, xlimit + 1, 1)], fontsize=12)
            ax.set_xlim(2 - 0.52, xlimit + 0.52)

        with mpl.rc_context({
            "legend.fontsize": 12, "ytick.labelsize": 12
        }):
            sns.histplot(
                df_truth_set,
                x="NumRepeatsAwayFromReference" if i == 0 else "MotifSize",
                hue=color_by,
                hue_order=hue_order,
                binwidth=1,
                multiple="stack" if not color_by else "fill",
                stat="proportion",
                discrete=True,
                legend=i == 1,
                palette=palette,
                ax=ax)

        if i == 0:
            ax.set_ylabel("Fraction of Alleles", fontsize=14)
        else:
            ax.set_ylabel("")

        if ax.get_legend():
            if color_by:
                sns.move_legend(ax, loc="upper right")
            if color_by == "Multiallelic":
                ax.get_legend().set_title("Multiallelic")
            elif color_by == "OverlapsSegDupIntervals":
                ax.get_legend().set_title("\n".join([
                    "TR Locus Overlaps ",
                    "Segmental Duplication",
                ]))
            ax.get_legend().get_title().set_horizontalalignment('center')

    if color_by:
        output_image_name += f".color_by_{color_by.lower().replace(' ', '_')}"

    output_path = os.path.join(args.output_dir, output_image_name + f".{args.image_type}")
    plt.savefig(output_path, bbox_extra_artists=extra_artists, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {output_path}")

    print(f"Plotted {len(df_truth_set):,d} allele records")


def plot_gene_info(df, args, excluding_introns_and_intergenic=False, use_MANE_genes=False):

    if use_MANE_genes:
        gene_region_column = "GeneRegionFromMane_V1"
        output_image_filename_suffix = ".MANE_v1"
        title_string = "MANE v1"
    else:
        gene_region_column = "GeneRegionFromGencode_V42"
        output_image_filename_suffix = ".gencode_v42"
        title_string = "Gencode v42"

    df.loc[:, gene_region_column] = df[gene_region_column].str.replace("CDS", "coding region")
    df.loc[:, gene_region_column] = df[gene_region_column].str.replace("exon", "exon of non-coding gene")

    if excluding_introns_and_intergenic:
        df = df[~df[gene_region_column].isin(["intron", "intergenic", "exon of non-coding gene"])]

    hue_order = ["intron", "exon of non-coding gene", "promoter", "5' UTR", "3' UTR", "coding region", "intergenic"]
    palette = sns.color_palette("tab20", n_colors=len(hue_order))
    palette = palette[1:] + palette[:1]

    if excluding_introns_and_intergenic:
        hue_order = hue_order[2:-1]
        palette = palette[2:-1]

    plt.rcParams.update({
        "legend.fontsize": 12,
        "ytick.labelsize": 12,
    })
    if not excluding_introns_and_intergenic:
        plt.rcParams.update({"legend.loc": "lower left" if use_MANE_genes else "upper left"})

    fig, ax = plt.subplots(1, 1, figsize=(8, 6) if args.width is None or args.height is None else (args.width, args.height))
    if args.show_title:
        suptitle_artist = fig.suptitle(f"Truth Set Alleles: Gene Region Overlap ({title_string})", fontsize=15)
        extra_artists = [suptitle_artist]
    else:
        extra_artists = []
    ax.xaxis.labelpad = ax.yaxis.labelpad = 15
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    xlimit = 30
    ax.set_xlabel("Motif Size (bp)", fontsize=13)
    ax.set_xticks(range(2, xlimit + 1, 1))
    ax.set_xticklabels([f"{x}" if x <= 6 or x % 3 == 0 else "" for x in range(2, xlimit + 1, 1)], fontsize=12)
    ax.set_xlim(2 - 0.52, xlimit + 0.52)

    sns.histplot(
        df,
        x="MotifSize",
        hue=gene_region_column,
        hue_order=hue_order,
        palette=palette,
        binwidth=1,
        multiple="fill",
        stat="proportion",
        discrete=True,
        legend=True,
        ax=ax)

    ax.set_ylabel("Fraction of Alleles", fontsize=14)

    if ax.get_legend():
        sns.move_legend(ax, loc=(1.05, 0.775 if excluding_introns_and_intergenic else 0.62))
        ax.get_legend().set_title("")
        ax.get_legend().set_frame_on(False)
        extra_artists.append(ax.get_legend())

    output_image_name = "truth_set_gene_overlap"
    if args.only_pure_repeats:
        output_image_name += ".only_pure_repeats"
    if excluding_introns_and_intergenic:
        output_image_name += ".excluding_introns_and_intergenic"
    else:
        output_image_name += ".all_regions"
    output_image_name += f"{output_image_filename_suffix}.{args.image_type}"

    output_path = os.path.join(args.output_dir, output_image_name)
    plt.savefig(output_path, bbox_extra_artists=extra_artists, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {output_path}")

    print(f"Plotted {len(df):,d} allele records")


def plot_allele_size_distribution_x3(df, args, color_by=None, hue_order=None, palette=None):

    """Plot allele size distribution split into 3 plots by motif size ranges"""
    fig, axes = plt.subplots(ncols=3, figsize=(30, 10) if args.width is None or args.height is None else (args.width, args.height), sharey=True)

    if args.show_title:
        fig.suptitle("Allele Size Distributions", fontsize=24)

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
        else:
            ax.set_ylabel(" ")

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
        ax.yaxis.grid(True, zorder=0, color="#D5D5F5")
        ax.yaxis.set_tick_params(labelleft=True)

        p = sns.histplot(
            df_current,
            x="NumBasePairsAwayFromReference",
            hue=color_by,
            hue_order=hue_order,
            palette=palette,
            multiple="stack",
            bins=[b-1.5 for b in range(-xlimit, xlimit + step_size, step_size)],
            zorder=3,
            ax=ax)

        p.set(yscale='log')
        p.set_title(f"{len(df_current):,d} alleles at {len(set(df_current.LocusId)):,d} loci ({label})", fontsize=18)
        if color_by:
            sns.move_legend(ax, loc="upper left")

        output_image_name = "allele_size_distribution_by_number_of_repeats.x3"
        if args.only_pure_repeats:
            output_image_name += ".only_pure_repeats"

        if color_by:
            output_image_name += f".color_by_{color_by.lower().replace(' ', '_')}"

        output_path = os.path.join(args.output_dir, output_image_name + f".{args.image_type}")
        fig.savefig(output_path, bbox_inches="tight", dpi=300)
        print(f"Saved {output_path}")

        print(f"Plotted {len(df):,d} allele records")


def plot_motif_distribution(df, args):

    print("Plotting allele distribution by motif size")

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
        f"-{hue_limit}", f"-{hue_limit} or more")
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].replace(
        f"+{hue_limit}", f"+{hue_limit} or more")

    plt.rcParams.update({
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
    })
    _, axes = plt.subplots(2, 1, figsize=(8, 12) if args.width is None or args.height is None else (args.width, args.height))

    title = "# of Repeats in hg38 at Truth Set Loci"

    ax = axes[0]
    p = sns.histplot(
        df,
        x="NumRepeatsInReference",
        binwidth=1,
        discrete=True,
        color="#1F77B4",
        zorder=3,
        ax=ax)

    ax.yaxis.grid(True, zorder=0, color="#D5D5F5")

    if args.show_title:
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

    xlimit = 50
    for ax in axes:
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xticks(np.arange(0, xlimit + 0.5, 5))
        ax.set_xlim(-0.5, xlimit + 0.5)

    ax.set_xlabel("# of Repeats in hg38", fontsize=14)
    ax.set_ylabel("Fraction of Alleles", fontsize=14)

    with mpl.rc_context({
        "text.usetex": True, "legend.fontsize": 12,
    }):
        ax.get_legend().set_title("\n".join([
            "\# of Repeats in",
            "Truth Set Allele",
            "Relative to hg38",
            "$_{contractions\ are\ <\ 0}$",
            "$_{expansions\ are\ >\ 0}$",
            "",
        ]), prop={'size': 14})
        sns.move_legend(ax, loc=(1.05, 0.15))
        ax.get_legend()._legend_box.align = "left"
        ax.get_legend().set_frame_on(False)

    output_image_name = "reference_locus_size_distribution"
    if args.only_pure_repeats:
        output_image_name += ".only_pure_repeats"

    output_path = os.path.join(args.output_dir, output_image_name + f".{args.image_type}")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    plt.close()
    print(f"Saved {output_path}")

    print(f"Plotted {len(df):,d} allele records")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--only-plot", choices=["1", "2", "3", "4", "5", "6", "7"], help="If specified, only plot this one figure")
    p.add_argument("--only-pure-repeats", action="store_true", help="If specified, only plot pure repeat alleles")
    p.add_argument("--width", type=float, help="Width of image")
    p.add_argument("--height", type=float, help="Height of image")
    p.add_argument("--show-title", action="store_true", help="If specified, show the title on the plots")
    p.add_argument("--image-type", choices=["svg", "png"], default="svg")

    p.add_argument("--output-dir", default=".")
    p.add_argument("--truth-set-alleles-table", default="../STR_truth_set.v1.alleles.tsv.gz")
    args = p.parse_args()

    input_table_path = args.truth_set_alleles_table

    print(f"Reading {input_table_path}")
    df = pd.read_table(input_table_path)

    print(f"Parsed {len(df):,d} rows from {input_table_path}")

    if args.only_pure_repeats:
        df = df[df.IsPureRepeat]
        print(f"Filtered to {len(df):,d} pure repeat alleles")

    df = df.rename(columns={"IsMultiallelic": "Multiallelic"})
    df["NumRepeatsAwayFromReference"] = df.NumRepeats - df.NumRepeatsInReference
    df["NumBasePairsAwayFromReference"] = df.NumRepeatsAwayFromReference * df.MotifSize
    df["Multiallelic"] = df["Multiallelic"].replace({True: "Yes", False: "No"})
    df["OverlapsSegDupIntervals"] = df["OverlapsSegDupIntervals"].replace({True: "Yes", False: "No"})

    if not args.only_plot or args.only_plot == "1":
        plot_allele_size_distribution(df, args, plot_type=1)
        plot_allele_size_distribution(df, args, plot_type=2)

    if not args.only_plot or args.only_plot == "2":
        df["Interruptions"] = df["IsPureRepeat"].replace({True: "No", False: "Yes"})
        #plot_allele_size_distribution(df, args, plot_type=1, color_by="Multiallelic", hue_order=["No", "Yes"], yaxis_fraction=True)
        #plot_allele_size_distribution(df, args, plot_type=3, color_by="Multiallelic", hue_order=["No", "Yes"])

        plot_allele_size_distribution(df, args, plot_type=1, color_by="Interruptions", hue_order=["No", "Yes"],
                                      yaxis_fraction=False, palette=["#1f77b4", "#8FD7E8"])
        #plot_allele_size_distribution(df, args, plot_type=3, color_by="Has Interruptions", hue_order=["No", "Yes"])

        #plot_allele_size_distribution(df, args, plot_type=1, color_by="OverlapsSegDupIntervals", hue_order=["No", "Yes"], yaxis_fraction=True)

    if not args.only_plot or args.only_plot == "3":
        plot_allele_size_distribution_x3(df, args)

    if not args.only_plot or args.only_plot == "4":
        plot_allele_size_and_motif_distribution(df, args, color_by="Multiallelic", hue_order=["No", "Yes"], palette=["#1F77B4", "#FF8C1A"])
        plot_allele_size_and_motif_distribution(df, args, color_by="OverlapsSegDupIntervals", hue_order=["No", "Yes"], palette=["#1F77B4", "#FF8C1A"])

    if not args.only_plot or args.only_plot == "5":
        plot_gene_info(df, args, excluding_introns_and_intergenic=False, use_MANE_genes=False)
        plot_gene_info(df, args, excluding_introns_and_intergenic=False, use_MANE_genes=True)
        plot_gene_info(df, args, excluding_introns_and_intergenic=True, use_MANE_genes=False)
        plot_gene_info(df, args, excluding_introns_and_intergenic=True, use_MANE_genes=True)

    if not args.only_plot or args.only_plot == "6":
        plot_motif_distribution(df, args)


if __name__ == "__main__":
    main()