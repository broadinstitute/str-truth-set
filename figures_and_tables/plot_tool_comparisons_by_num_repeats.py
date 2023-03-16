import argparse
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import pandas as pd
import seaborn as sns
import sys


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})

GREEN_COLOR = "#50AA44"

NO_CALL_LABEL = "No Call"
FILTERED_CALL_LABEL = "Filtered"
HOM_REF_LABEL = "Called Hom Ref"
HET_REF_LABEL = "Called Het Ref"
WRONG_DIRECTION_LABEL = "Wrong Direction"


def bin_num_repeats_wrapper(bin_size=1):
    """Create a function that converts the integer number of STR repeats to a histogram bin"""

    def bin_num_repeats(row):
        num_repeats = row["DiffFromRefRepeats: Allele: Truth"]
        if num_repeats == 0:
            return "0"

        sign = "-" if num_repeats < 0 else "+"
        num_repeats = abs(num_repeats)
        if row["Coverage"] == "exome":
            if 13 <= num_repeats:
                return f"{sign}13 or more"
        else:
            if 21 <= num_repeats <= 25:
                return f"{sign}21 to {sign}25"
            if 26 <= num_repeats <= 30:
                return f"{sign}26 to {sign}30"
            if 31 <= num_repeats:
                return f"{sign}31 or more"

        if bin_size == 1:
            return f"{sign}{int(num_repeats)}"

        start = int((num_repeats - 1)/bin_size) * bin_size
        end = start + bin_size
        return f"{sign}{int(start + 1)} to {sign}{int(end)}"

    return bin_num_repeats


def bin_repeats_vs_truth_wrapper():
    """Create a function that converts the difference between a tool's output allele size (number of STR repeats) and
    the true number of repeats to a histogram bin"""

    def bin_repeats_vs_truth(num_repeats_diff):
        if pd.isna(num_repeats_diff):
            return NO_CALL_LABEL

        sign = "-" if num_repeats_diff < 0 else ""
        num_repeats_diff = abs(num_repeats_diff)
        if num_repeats_diff <= 2:
            return f"{sign}{int(num_repeats_diff)}"
        elif 3 <= num_repeats_diff <= 7:
            return f"{sign}3 to {sign}7"
        elif 8 <= num_repeats_diff <= 20:
            return f"{sign}8 to {sign}20"
        elif num_repeats_diff >= 21:
            return f"{sign}21" + " or more"
        else:
            raise ValueError(f"Unexpected num_repeats_diff: {num_repeats_diff}")

    return bin_repeats_vs_truth


def define_hue_column(df, tool):
    """Adds a new column to the input table which is used to color the plots"""
    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = df[f"DiffRepeats: Allele: {tool} - Truth"].apply(
        bin_repeats_vs_truth_wrapper())

    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        (df["DiffFromRefRepeats: Allele: Truth"] != 0) &
        ~df[f"DiffFromRefRepeats: Allele: Truth"].isna() &
        ~df[f"DiffFromRefRepeats: Allele: {tool}"].isna() & (
        df[f"DiffFromRefRepeats: Allele: {tool}"] != 0) & (
            np.sign(df[f"DiffFromRefRepeats: Allele: Truth"]) != np.sign(df[f"DiffFromRefRepeats: Allele: {tool}"])
        ),
        WRONG_DIRECTION_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        (df["DiffFromRefRepeats: Allele: Truth"] != 0) & ~df[f"IsRef: Allele: {tool}"].isna() & df[f"IsRef: Allele: {tool}"],
        HET_REF_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        (df["DiffFromRefRepeats: Allele: Truth"] != 0) & ~df[f"IsHomRef: {tool}"].isna() & df[f"IsHomRef: {tool}"],
        HOM_REF_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])


DISTRIBUTION_BY_NUM_REPEATS_FIGURE_SIZE = (20, 9)


def plot_empty_image(figure_title, message):
    fig, ax = plt.subplots(1, 1, figsize=DISTRIBUTION_BY_NUM_REPEATS_FIGURE_SIZE, dpi=80)
    fig.suptitle(figure_title, fontsize=20)
    ax.axis('off')
    text = ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=24)

    plt.gcf().canvas.draw()

    bbox = text.get_window_extent()
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    rect = patches.Rectangle(bbox.min, bbox.width, bbox.height, linewidth=1, edgecolor='black', facecolor='none', ls='dashed')
    ax.add_patch(rect)


def plot_distribution_by_num_repeats(
    df,
    x_column,
    hue_column,
    hue_order=None,
    tool_name=None,
    palette=None,
    figure_title=None,
):

    fig, axes = plt.subplots(1, 2, figsize=DISTRIBUTION_BY_NUM_REPEATS_FIGURE_SIZE, dpi=80)

    if figure_title:
        fig.suptitle(figure_title, fontsize=20)

    for i, ax in enumerate(axes):
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15

        ax.set_xlabel("True Allele Size Minus Number of Repeats in Reference Genome", fontsize=16)
        if i == 0:
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.set_ylabel("Number of Alleles", fontsize=16)
        else:
            ax.set_ylabel("Fraction of Alleles", fontsize=16)

        ax.tick_params(axis="x")

        sns.histplot(
            df,
            x=x_column,
            hue=hue_column,
            hue_order=hue_order,
            binwidth=1,
            palette=palette,
            multiple="stack" if i == 0 else "fill",
            stat="count" if i == 0 else "proportion",
            discrete=True,
            legend=i == 0,
            ax=axes[i])

        fig.tight_layout()
        axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45, horizontalalignment="right", rotation_mode="anchor", fontsize=11)
        axes[i].set_yticklabels(axes[i].get_yticklabels(), fontsize=12)
        if i == 1:
            # add n=.. above each bar
            n_lookup = dict(df.groupby(x_column).count().LocusId)
            for j, (xtick, text) in enumerate(zip(axes[i].get_xticks(), axes[i].get_xticklabels())):
                ax.text(xtick, 1.018, f"{n_lookup[text.get_text()]:,d} alleles", ha="left", va="bottom", color="#777777", rotation=45)

    if tool_name:
        l = axes[0].get_legend()
        l.set_title(f"{tool_name} Call\nvs\nTrue Allele Size\n")
        axes[0].get_legend().get_title().set_horizontalalignment('center')
        axes[0].get_legend().set_frame_on(False)
    else:
        axes[0].get_legend().set_title(f"")

    fig.tight_layout()


def hue_sorter(value):
    """Defines the sort order of colors in the plot_distribution_by_num_repeats plot"""
    if value == FILTERED_CALL_LABEL:
        return -2500
    elif value == NO_CALL_LABEL:
        return -2000
    elif value == HOM_REF_LABEL:
        return -1500
    elif value == HET_REF_LABEL:
        return -1000
    elif value == WRONG_DIRECTION_LABEL:
        return -500
    else:
        return int(float(value.split(" ")[0]))


def generate_all_distribution_by_num_repeats_plots(df, output_image_dir, max_plots=None):
    plot_counter = 0
    for motif_size in ("all_motifs", "STR", "TR", "TR25+"):
        for pure_repeats in (True, False,):
            for coverage in ("40x", "30x", "20x", "10x", "05x", "exome", ):
                for tool_label in ("ExpansionHunter", "GangSTR", "HipSTR"):  #"GangSTR__Q_over_0.8", "ExpansionHunter_Filtered":
                    for genotype_subset in ("all_genotypes", "HET", "HOM", "MULTI"):
                        if max_plots is not None and plot_counter >= max_plots:
                            print(f"Exiting after generating {plot_counter} plot(s)")
                            sys.exit(0)

                        print("-"*100)

                        figure_title_line1 = ""
                        figure_title_line2 = ""

                        df2 = df.copy()
                        df2 = df2[(df2["Coverage"] == coverage)]

                        if coverage == "exome":
                            df2 = df2[~df2["Genotype: GangSTR"].isna() | ~df2["Genotype: ExpansionHunter"].isna()]
                            df2 = df2[~df2["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

                        filter_description = []
                        output_image_filename = "tool_accuracy_by_true_allele_size"
                        if motif_size == "all_motifs":
                            filter_description.append("all motif sizes")
                            output_image_filename += ".all_motifs"
                            df2 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
                        elif motif_size == "STR":
                            filter_description.append("2bp to 6bp motifs")
                            output_image_filename += ".2to6bp_motifs"
                            df2 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 6)]
                        elif motif_size == "TR":
                            filter_description.append(f"7bp to 24bp motifs")
                            output_image_filename += ".7to24bp_motifs"
                            df2 = df2[(7 <= df2["MotifSize"]) & (df2["MotifSize"] <= 24)]
                        elif motif_size == "TR25+":
                            filter_description.append(f"25bp to 50bp motifs")
                            output_image_filename += ".25to50bp_motifs"
                            df2 = df2[(25 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
                        else:
                            raise ValueError(f"Unexpected motif_size value: {motif_size}")

                        if genotype_subset == "all_genotypes":
                            output_image_filename += ".all_genotypes"
                        elif genotype_subset == "HET":
                            filter_description.append("HET")
                            output_image_filename += ".HET"
                            df2 = df2[df2["SummaryString"].str.contains(":HET")]
                        elif genotype_subset == "HOM":
                            filter_description.append(f"HOM")
                            output_image_filename += ".HOM"
                            df2 = df2[df2["SummaryString"].str.contains(":HOM")]
                        elif genotype_subset == "MULTI":
                            filter_description.append(f"multi-allelic")
                            output_image_filename += ".MULTI"
                            df2 = df2[df2["SummaryString"].str.contains(":MULTI")]
                        else:
                            raise ValueError(f"Unexpected motif_size value: {motif_size}")

                        if pure_repeats:
                            filter_description.append("uninterrupted repeats")
                            output_image_filename += ".pure_repeats"
                            df2 = df2[df2["IsPureRepeat"]]
                        else:
                            filter_description.append("repeats with interruptions")
                            output_image_filename += ".with_interruptions"
                            df2 = df2[~df2["IsPureRepeat"]]

                        output_image_filename += f".{coverage}"

                        if "GangSTR__Q_over_0.8" in tool_label:
                            tool = "GangSTR"
                            filter_description.append(f"filtered to Q > 0.8")
                            output_image_filename += f".{tool}_filtered"
                        elif "ExpansionHunter_Filtered" in tool_label:
                            tool = "ExpansionHunter"
                            filter_description.append(f"filtered")
                            output_image_filename += f".{tool}_filtered"
                        else:
                            tool = tool_label
                            output_image_filename += f".{tool}"

                        coverage_label = f"exome data" if coverage == "exome" else f"{coverage} genome data"

                        # skip_condition1: HipSTR doesn't support motifs larger than 9bp
                        skip_condition1 = tool == "HipSTR" and motif_size != "STR"
                        # skip_condition2: not enough alleles to draw a histogram
                        skip_condition2 = len(df2) < 10

                        if skip_condition1 or skip_condition2:
                            if skip_condition1:
                                message = "HipSTR doesn't support motif sizes larger than 9bp"
                            elif skip_condition2:
                                message = "Not enough alleles to create plot"

                            figure_title_line1 += f"{tool} {coverage_label}"
                            print(figure_title_line1)
                            #if len(df2) > 0:
                            #    figure_title_line2 += f"{len(df2):,d} alleles at "
                            figure_title_line2 += f"{len(set(df2.LocusId)):,d} loci ("+", ".join(filter_description)+")"

                            print(f"Skipping..  {message}")
                            plot_empty_image(figure_title_line1 + "\n\n" + figure_title_line2, message)
                            plt.savefig(f"{output_image_dir}/{output_image_filename}.svg")
                            print(f"Saved {output_image_dir}/{output_image_filename}.svg")
                            plt.close()
                            continue

                        num_alleles_exactly_right = sum(df2[f"DiffRepeats: Allele: {tool} - Truth"] == 0)
                        print("Keeping", len(df2), "out of ", len(df), "rows")
                        figure_title_line1 += f"{tool} got {num_alleles_exactly_right:,d} out of {len(df2):,d} alleles ({100*num_alleles_exactly_right/len(df2):0.1f}%) exactly right for {coverage_label}"
                        print(figure_title_line1)
                        figure_title_line2 += f"Showing results for {len(set(df2.LocusId)):,d} loci (" + ", ".join(filter_description) + f")"
                        print(figure_title_line2)

                        print(tool, sum(df2[f"DiffFromRefRepeats: Allele: {tool}"].isna()), "out of", len(df2), f"{tool} - Ref' values are NaN")
                        print(tool, sum(df2[f"DiffRepeats: Allele: {tool} - Truth"].isna()), "out of", len(df2), f"{tool} - Truth' values are NaN")

                        if "GangSTR__Q_over_0.8" in tool_label:
                            df2.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
                                ~df2[f"Q: GangSTR"].isna() & (df2[f"Q: GangSTR"].astype(float) <= 0.8),
                                "Filtered",
                                df2[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

                        elif "ExpansionHunter_Filtered" in tool_label:
                            df2.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
                                ~df2[f"CI size: Allele: ExpansionHunter"].isna() & (
                                        df2[f"CI size: Allele: ExpansionHunter"].astype(float)/df2["NumRepeats: Allele: ExpansionHunter"] > 2),
                                "Filtered",
                                df2[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

                        hue_values = set(df2.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"])
                        print("Hue values:", sorted(hue_values, key=hue_sorter))

                        palette = []
                        if FILTERED_CALL_LABEL in hue_values:
                            palette.append("#f5f5f5")
                        if NO_CALL_LABEL in hue_values:
                            palette.append("#c9c9c9")
                        if HOM_REF_LABEL in hue_values:
                            palette.append("#C9342A")
                        if HET_REF_LABEL in hue_values:
                            palette.append("#734B4B")
                        if WRONG_DIRECTION_LABEL in hue_values:
                            palette.append("#99FFEF")

                        n_blue_colors = sum(1 for h in hue_values if h.startswith("-"))
                        n_orange_colors = len(hue_values) - n_blue_colors - len(palette) - 1

                        try:
                            palette += list(sns.color_palette("Blues_r", n_colors=n_blue_colors))
                            palette += [GREEN_COLOR]
                            palette += list(sns.color_palette("Oranges_r", n_colors=n_orange_colors)[::-1])

                            plot_distribution_by_num_repeats(
                                df2,
                                x_column="DiffFromRefRepeats: Allele: Truth (bin)",
                                hue_column=f"DiffRepeats: Allele: {tool} - Truth (bin)",
                                hue_order=sorted(hue_values, key=hue_sorter),
                                tool_name=tool,
                                palette=palette,
                                figure_title=figure_title_line1 + "\n\n" + figure_title_line2
                            )

                            for ext in ".svg", ".png":
                                plt.savefig(f"{output_image_dir}/{output_image_filename}{ext}")
                                print(f"Saved plot #{plot_counter+1}: {output_image_dir}/{output_image_filename}{ext}")
                            plt.close()

                        except Exception as e:
                            print(f"ERROR: {e}")
                            #import traceback
                            #traceback.print_exc()

                        plot_counter += 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="If specified, only generate this many plots. Useful for testing")
    p.add_argument("--output-dir", default="../tool_comparison/figures")
    p.add_argument("--verbose", action="store_true", help="Print additional info")
    p.add_argument("combined_tool_results_tsv", nargs="?", default="../tool_comparison/combined.results.alleles.tsv.gz")
    args = p.parse_args()

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)
    df = df[df["IsFoundInReference"] & (df["PositiveOrNegative"] == "positive")]

    if args.verbose:
        print("Num loci:")
        print(df.groupby(["PositiveOrNegative", "Coverage", "IsPureRepeat"]).count().LocusId/2)

    print("Computing additional columns...")
    df.loc[:, "DiffFromRefRepeats: Allele: Truth (bin)"] = df.apply(bin_num_repeats_wrapper(bin_size=2), axis=1)
    for tool in ("ExpansionHunter", "GangSTR", "HipSTR",):
        define_hue_column(df, tool)

    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    print("Generating plots...")
    generate_all_distribution_by_num_repeats_plots(df, args.output_dir, max_plots=args.n)

    print(f"Done")


if __name__ == "__main__":
    main()