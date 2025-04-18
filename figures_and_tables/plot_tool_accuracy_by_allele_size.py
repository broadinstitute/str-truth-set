import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

from matplotlib import patches

sns.set_context(font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
})

GREEN_COLOR = "#50AA44"

NO_CALL_LABEL = "No Call"
FILTERED_CALL_LABEL = "Filtered"
HOM_REF_LABEL = "Called Hom Ref"
HET_REF_LABEL = "Called Het Ref"
WRONG_DIRECTION_LABEL = "Wrong Direction"
NO_DIFFERENCE_LABEL = "Same"

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


def bin_repeats_vs_truth_wrapper(tool):
    """Create a function that converts the difference between a tool's output allele size (number of STR repeats) and
    the true number of repeats to a histogram bin"""

    def bin_repeats_vs_truth(row):
        num_repeats_diff = row[f"DiffRepeats: Allele: {tool} - Truth"]
        if pd.isna(num_repeats_diff):
            return NO_CALL_LABEL

        num_repeats = int(row[f"NumRepeats: Allele: {tool}"])
        sign = "-" if num_repeats_diff < -1 else ""
        num_repeats_diff = abs(num_repeats_diff)
        if (num_repeats_diff <= 1) or (num_repeats > 120 and num_repeats_diff <= 2) or (num_repeats > 240 and num_repeats_diff <= 3) or (num_repeats > 360 and num_repeats_diff <= 4):
            return NO_DIFFERENCE_LABEL
        elif num_repeats_diff == 2:
            return f"{sign}{int(num_repeats_diff)}"
        elif 3 <= num_repeats_diff <= 4:
            return f"{sign}3 to {sign}4"
        elif 5 <= num_repeats_diff <= 7:
            return f"{sign}5 to {sign}7"
        elif 8 <= num_repeats_diff <= 20:
            return f"{sign}8 to {sign}20"
        elif num_repeats_diff >= 21:
            return f"{sign}21" + " or more"
        else:
            raise ValueError(f"Unexpected num_repeats_diff: {num_repeats_diff}")

    return bin_repeats_vs_truth


def define_hue_column(df, tool):
    """Adds a new column to the input table which is used to color the plots"""
    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = df.apply(bin_repeats_vs_truth_wrapper(tool), axis=1)

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


def plot_empty_image(figure_title, message, width, height):
    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    fig.suptitle(figure_title, fontsize=17)
    ax.axis('off')
    text = ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=17)

    plt.gcf().canvas.draw()

    bbox = text.get_window_extent()
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    rect = patches.Rectangle(bbox.min, bbox.width, bbox.height, linewidth=1, edgecolor='black', facecolor='none', ls='dashed')
    ax.add_patch(rect)


def plot_accuracy_by_allele_size(
    df,
    args,
    x_column,
    hue_column,
    hue_order=None,
    tool_name=None,
    palette=None,
    figure_title=None,
):

    fig, axes = plt.subplots(1, 2, figsize=(args.width, args.height))

    if figure_title:
        fig.suptitle(figure_title, fontsize=17)

    for i, ax in enumerate(axes):
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15

        ax.set_xlabel("True Allele Size Minus Number of Repeats in Reference Genome", fontsize=14)
        if i == 0:
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.set_ylabel("Number of Alleles", fontsize=14)
        else:
            ax.set_ylabel("Fraction of Alleles", fontsize=14)

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
                ax.text(xtick, 1.018, f"{n_lookup[text.get_text()]:,d}", ha="left", va="bottom", color="#777777", fontsize=12, rotation=45)
            ax.text(0, 1.15, "Alleles Per Bin", ha="left", va="bottom", color="#777777", fontsize=12, transform=ax.transAxes)

    if tool_name:
        l = axes[0].get_legend()
        l.set_title(f"{tool_name} Call\nvs\nTrue Allele Size\n")
        axes[0].get_legend().get_title().set_horizontalalignment('center')
        axes[0].get_legend().set_frame_on(False)
    else:
        axes[0].get_legend().set_title(f"")

    fig.tight_layout()


def hue_sorter(value):
    """Defines the sort order of colors in the plot_distribution_by_allele_size plot"""
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
    elif value == NO_DIFFERENCE_LABEL:
        return 0
    else:
        return int(float(value.split(" ")[0]))


def generate_all_plots(df, args):
    arg_string = " ".join(sys.argv)
    plot_counter = 0

    start_with_plot_i = args.start_with_plot_i
    if start_with_plot_i is None or start_with_plot_i < 0:
        start_with_plot_i = 0

    max_plots = args.n
    if max_plots is None:
        max_plots = 10**9

    for tool in ["ExpansionHunter", "ExpansionHunter-dev", "GangSTR", "HipSTR", "constrain"] if not args.tool else [args.tool]:
        for q_threshold in [0, 0.1, 0.5, 0.9] if args.q_threshold is None else [args.q_threshold]:
            if q_threshold > 0 and not args.only_print_total_number_of_plots and tool != "constrain":
                df2 = df.copy()
                q_column = f"Q: Allele: {tool}" if tool in ("ExpansionHunter", "ExpansionHunter-dev") else f"Q: {tool}"
                df2.loc[df2[q_column] < q_threshold, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = FILTERED_CALL_LABEL
            else:
                df2 = df

            for coverage in ["40x", "30x", "20x", "10x", "exome"] if not args.coverage else [args.coverage]:
                if not args.only_print_total_number_of_plots:
                    df3 = df2[(df2["Coverage"] == coverage)]
                else:
                    df3 = df2

                if coverage == "exome" and not args.only_print_total_number_of_plots:
                    df3 = df3[~df3["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

                for motif_size in ["all_motifs", "2bp", "3-6bp", "7-24bp", "25+bp"] if args.min_motif_size is None and args.max_motif_size is None else [f"{args.min_motif_size}-{args.max_motif_size}bp"]:
                    if not args.only_print_total_number_of_plots:
                        if motif_size == "all_motifs":
                            df4 = df3[(2 <= df3["MotifSize"]) & (df3["MotifSize"] <= 50)]
                        elif motif_size == "2bp":
                            df4 = df3[df3["MotifSize"] == 2]
                        elif motif_size == "3-6bp":
                            df4 = df3[(3 <= df3["MotifSize"]) & (df3["MotifSize"] <= 6)]
                        elif motif_size == "7-24bp":
                            df4 = df3[(7 <= df3["MotifSize"]) & (df3["MotifSize"] <= 24)]
                        elif motif_size == "25+bp":
                            df4 = df3[(25 <= df3["MotifSize"]) & (df3["MotifSize"] <= 50)]
                        else:
                            if args.min_motif_size is None and args.max_motif_size is None:
                                raise ValueError(f"Unexpected motif_size value: {motif_size}")
                            df4 = df3
                            if args.min_motif_size:
                                df4 = df4[(args.min_motif_size <= df4["MotifSize"])]
                            if args.max_motif_size:
                                df4 = df4[(df4["MotifSize"] <= args.max_motif_size)]

                    else:
                        df4 = df3

                    for genotype_subset in ["all", "HET", "HOM", "MULTI"] if not args.genotype else [args.genotype]:
                        if not args.only_print_total_number_of_plots:
                            if genotype_subset == "all":
                                df5 = df4
                            elif genotype_subset == "HET":
                                df5 = df4[df4["SummaryString"].str.contains(":HET")]
                            elif genotype_subset == "HOM":
                                df5 = df4[df4["SummaryString"].str.contains(":HOM")]
                            elif genotype_subset == "MULTI":
                                df5 = df4[df4["SummaryString"].str.contains(":MULTI")]
                            else:
                                raise ValueError(f"Unexpected genotype_subset value: {genotype_subset}")
                        else:
                            df5 = df4

                        for pure_repeats in ["both", True, False] if "--only-pure-repeats" not in arg_string else [True]:
                            if pure_repeats == "both" or args.only_print_total_number_of_plots:
                                df6 = df5
                            elif pure_repeats:
                                df6 = df5[df5["IsPureRepeat"]]
                            else:
                                df6 = df5[~df5["IsPureRepeat"]]

                            hide_no_call_loci_option = "--hide-no-call-loci" in arg_string
                            show_no_call_loci_option = "--show-no-call-loci" in arg_string
                            for only_loci_with_calls_by_this_tool in [False, True] if not hide_no_call_loci_option and not show_no_call_loci_option else ([True] if args.hide_no_call_loci else [False]):
                                if only_loci_with_calls_by_this_tool and not args.only_print_total_number_of_plots:
                                    df7 = df6[df6[f"DiffRepeats: Allele: {tool} - Truth (bin)"] != NO_CALL_LABEL]
                                else:
                                    df7 = df6

                                for exclude_filtered_loci in [False, True] if q_threshold > 0 else [False]:
                                    if exclude_filtered_loci and not args.only_print_total_number_of_plots:
                                        df_plot = df7[df7[f"DiffRepeats: Allele: {tool} - Truth (bin)"] != FILTERED_CALL_LABEL]
                                    else:
                                        df_plot = df7

                                    if plot_counter < start_with_plot_i or args.only_print_total_number_of_plots:
                                        plot_counter += 1
                                        continue

                                    print("-"*100)

                                    figure_title_line1 = ""
                                    figure_title_line2 = ""

                                    filter_description = []
                                    output_image_filename = "tool_accuracy_by_true_allele_size"
                                    if motif_size == "all_motifs":
                                        filter_description.append("all motif sizes")
                                        output_image_filename += ".all_motifs"
                                    elif motif_size == "2bp":
                                        filter_description.append("2bp motifs")
                                        output_image_filename += ".2bp_motifs"
                                    elif motif_size == "3-6bp":
                                        filter_description.append("3bp to 6bp motifs")
                                        output_image_filename += ".3to6bp_motifs"
                                    elif motif_size == "7-24bp":
                                        filter_description.append(f"7bp to 24bp motifs")
                                        output_image_filename += ".7to24bp_motifs"
                                    elif motif_size == "25+bp":
                                        filter_description.append(f"25bp to 50bp motifs")
                                        output_image_filename += ".25to50bp_motifs"
                                    else:
                                        if args.min_motif_size is None and args.max_motif_size is None:
                                            raise ValueError(f"Unexpected motif_size value: {motif_size}")
                                        filter_description.append(f"{args.min_motif_size}bp to {args.max_motif_size}bp motifs")
                                        output_image_filename += "." + motif_size.replace("-", "to") + "_motifs"

                                    if genotype_subset == "all":
                                        output_image_filename += ".all_genotypes"
                                    elif genotype_subset == "HET":
                                        filter_description.append("HET")
                                        output_image_filename += ".HET"
                                    elif genotype_subset == "HOM":
                                        filter_description.append(f"HOM")
                                        output_image_filename += ".HOM"
                                    elif genotype_subset == "MULTI":
                                        filter_description.append(f"multi-allelic")
                                        output_image_filename += ".MULTI"
                                    else:
                                        raise ValueError(f"Unexpected genotype_subset value: {genotype_subset}")

                                    if pure_repeats == "both":
                                        pass
                                    elif pure_repeats:
                                        filter_description.append("pure repeats")
                                        output_image_filename += ".pure_repeats"
                                    else:
                                        filter_description.append("repeats with interruptions")
                                        output_image_filename += ".with_interruptions"

                                    output_image_filename += f".{coverage}"

                                    if only_loci_with_calls_by_this_tool:
                                        filter_description.append(f"exclude no-call loci")
                                        output_image_filename += f".exclude_no_call_loci"

                                    if exclude_filtered_loci:
                                        filter_description.append(f"exclude filtered loci")
                                        output_image_filename += f".exclude_filtered_loci"

                                    if q_threshold > 0:
                                        filter_description.append(f"Q ≥ {q_threshold}")
                                        output_image_filename += f".Q{int(100*q_threshold)}"

                                    output_image_filename += f".{tool}"

                                    coverage_label = f"exome data" if coverage == "exome" else f"{coverage} genome data"

                                    # skip_condition1: HipSTR doesn't support motifs larger than 9bp
                                    skip_condition1 = tool == "HipSTR" and (motif_size not in ("2bp", "3-6bp")
                                            and (args.min_motif_size is None or args.min_motif_size > 9)
                                            and (args.max_motif_size is None or args.max_motif_size > 9)
                                    )
                                    # skip_condition2: not enough alleles to draw a histogram
                                    skip_condition2 = len(df_plot) < 10

                                    n_locus_ids = len(set(df_plot.LocusId))
                                    if skip_condition1 or skip_condition2:
                                        if skip_condition1:
                                            message = "HipSTR doesn't support motif sizes larger than 9bp"
                                        elif skip_condition2:
                                            message = "Not enough alleles to create plot"

                                        figure_title_line1 += f"{tool} {coverage_label}"
                                        figure_title_line2 += f"{n_locus_ids:,d} loci ("+", ".join(filter_description)+")"
                                        print(figure_title_line1)
                                        print(figure_title_line2)
                                        print(f"Skipping..  {message}")
                                        plot_empty_image(figure_title_line1 + "\n\n" + figure_title_line2, message, args.width, args.height)

                                        output_path = os.path.join(args.output_dir, f"{output_image_filename}.{args.image_type}")
                                        plt.savefig(output_path, dpi=300)
                                        print(f"Saved {output_path}")
                                        plt.close()
                                        continue

                                    print(f"Generating plot #{plot_counter}")
                                    num_alleles_exactly_right = sum(df_plot[f"DiffRepeats: Allele: {tool} - Truth (bin)"] == "0") + sum(df_plot[f"DiffRepeats: Allele: {tool} - Truth (bin)"] == NO_DIFFERENCE_LABEL)
                                    hue_values = set(df_plot.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"])

                                    figure_title_line1 += f"{tool} got {num_alleles_exactly_right:,d} out of {len(df_plot):,d} alleles ({100*num_alleles_exactly_right/len(df_plot):0.1f}%) exactly right for {coverage_label}"
                                    figure_title_line2 += f"Showing results for {n_locus_ids:,d} loci (" + ", ".join(filter_description) + f")"

                                    print(figure_title_line1)
                                    print(figure_title_line2)
                                    print(f"Plotting {len(df_plot):,d} out of {len(df):,d} rows")
                                    print("Hue values:", sorted(hue_values, key=hue_sorter))

                                    palette = []
                                    if FILTERED_CALL_LABEL in hue_values:
                                        palette.append("#f5f5f5")
                                    if NO_CALL_LABEL in hue_values:
                                        palette.append("#c5c5c5")
                                    if HOM_REF_LABEL in hue_values:
                                        palette.append("#C9342A")
                                    if HET_REF_LABEL in hue_values:
                                        palette.append("#734B4B")
                                    if WRONG_DIRECTION_LABEL in hue_values:
                                        palette.append("#99FFEF")

                                    n_blue_colors = sum(1 for h in hue_values if h.startswith("-"))
                                    if n_blue_colors > 0:
                                        palette += list(sns.color_palette("Blues_r", n_colors=n_blue_colors))
                                    if '0' in hue_values or NO_DIFFERENCE_LABEL in hue_values:
                                        palette += [GREEN_COLOR]
                                    n_orange_colors = len(hue_values) - len(palette)
                                    if n_orange_colors > 0:
                                        palette += list(sns.color_palette("Oranges_r", n_colors=n_orange_colors)[::-1])

                                    try:
                                        plot_accuracy_by_allele_size(
                                            df_plot,
                                            args,
                                            x_column="DiffFromRefRepeats: Allele: Truth (bin)",
                                            hue_column=f"DiffRepeats: Allele: {tool} - Truth (bin)",
                                            hue_order=sorted(hue_values, key=hue_sorter),
                                            tool_name=tool,
                                            palette=palette,
                                            figure_title=figure_title_line1 + "\n\n" + figure_title_line2 if args.show_title else None
                                        )

                                        output_path = os.path.join(args.output_dir, f"{output_image_filename}.{args.image_type}")
                                        plt.savefig(output_path, dpi=300)
                                        print(f"Saved plot #{plot_counter}: {output_path}")
                                        plt.close()

                                    except Exception as e:
                                        print(f"ERROR: {e}")

                                    plot_counter += 1
                                    if plot_counter >= start_with_plot_i + max_plots:
                                        print(f"Exiting after generating {plot_counter-start_with_plot_i} plot(s)")
                                        return

    if args.only_print_total_number_of_plots:
        print(f"Total: {plot_counter:,d} plots")



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start-with-plot-i", type=int, help="If specified, start with this plot number")
    p.add_argument("-n", type=int, help="If specified, only generate this many plots")
    p.add_argument("--width", default=20, type=float)
    p.add_argument("--height", default=9, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"], help="Image type to generate")
    p.add_argument("--show-title", action="store_true", help="Show plot title")

    p.add_argument("--output-dir", default=".")
    p.add_argument("--verbose", action="store_true", help="Print additional info")
    p.add_argument("--only-print-total-number-of-plots", action="store_true", help="Don't generate any plots. Only "
                   "print the total number of plots that would be generated.")

    g = p.add_argument_group("Filters")
    g.add_argument("--tool", choices={
        "ExpansionHunter", "ExpansionHunter-dev", "GangSTR", "HipSTR", "constrain", "TRGT", "LongTR", "NewTruthSet"},
        help="Plot only this tool")
    g.add_argument("--q-threshold", type=float, help="Plot only this Q threshold")
    g.add_argument("--coverage", help="Plot only this coverage (example: \"20x\" or \"exome\")")
    g.add_argument("--min-motif-size", type=int, help="Min motif size")
    g.add_argument("--max-motif-size", type=int, help="Max motif size")
    g.add_argument("--genotype", choices=["all", "HET", "HOM", "MULTI"], help="Plot only this genotype")
    g.add_argument("--only-pure-repeats", action="store_true", help="Plot only loci with pure repeats")
    g2 = g.add_mutually_exclusive_group()
    g2.add_argument("--show-no-call-loci", action="store_true", help="Show loci with no call")
    g2.add_argument("--hide-no-call-loci", action="store_true", help="Hide loci with no call")

    p.add_argument("combined_tool_results_tsv", nargs="?", default="../tool_comparison/combined.results.alleles.tsv.gz")
    args = p.parse_args()

    # print command line args from argparse
    print("Command line args:")
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")

    if args.only_print_total_number_of_plots:
        df = pd.DataFrame()
        generate_all_plots(df, args)
        return

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)

    before = len(df)
    df = df[~df["DiffFromRefRepeats: Allele: Truth"].isna()]
    print(f"Filtered out {before - len(df)} out of {before} ({(before - len(df))/before:0.1%}) rows with no NumRepeats")

    if "IsFoundInReference" in df.columns:
        df = df[df["IsFoundInReference"]]
    else:
        print("WARNING: IsFoundInReference column not found in input file. Assuming all loci are found in reference...")
        df["IsFoundInReference"] = True

    if "PositiveOrNegative" in df.columns:
        df = df[df["PositiveOrNegative"] == "positive"]
    else:
        print("WARNING: PositiveOrNegative column not found in input file. Assuming all loci are positive...")
        df["PositiveOrNegative"] = "positive"

    if "Coverage" not in df.columns:
        df["Coverage"] = args.coverage
    if f"Q: {args.tool}" not in df.columns:
        df[f"Q: {args.tool}"] = 1

    if args.verbose:
        print("Num loci:")
        print(df.groupby(["PositiveOrNegative", "Coverage", "IsPureRepeat"]).count().LocusId/2)

    print("Computing additional columns...")
    df.loc[:, "DiffFromRefRepeats: Allele: Truth (bin)"] = df.apply(bin_num_repeats_wrapper(bin_size=2), axis=1)
    for tool in ([args.tool] or ("ExpansionHunter", "ExpansionHunter-dev", "GangSTR", "HipSTR", "constrain", "TRGT", "LongTR")):
        define_hue_column(df, tool)

    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    print(f"Generating {str(args.n) + ' ' if args.n else ''}plots",
          f"starting with plot #{args.start_with_plot_i}" if args.start_with_plot_i else "")
    generate_all_plots(df, args)

    print(f"Done")


if __name__ == "__main__":
    main()
