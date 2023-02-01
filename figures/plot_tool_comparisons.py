import argparse
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
import pandas as pd
import seaborn as sns
import sys

sns.set_context("notebook", font_scale=1.1, rc={
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
        sign = "-" if num_repeats < 0 else "+"
        num_repeats = abs(num_repeats)

        if num_repeats == 0:
            return "0"

        if row["coverage"] == "exome":
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


def bin_repeat_size_bp_wrapper(threshold=55, bin_size=1):
    """Create a function that converts the integer STR repeat size in base pairs to a histogram bin"""
    threshold = threshold - threshold % bin_size + bin_size + 1
    def bin_repeat_size_bp(repeat_size):
        sign = "-" if repeat_size < 0 else ""
        repeat_size = abs(repeat_size)
        if repeat_size == 0:
            return "0"
        elif repeat_size >= threshold:
            return f"{sign}{threshold}bp or more"
        else:
            if bin_size > 1:
                start = int((repeat_size - 1)/bin_size) * bin_size
                end = start + bin_size
                return f"{sign}{start + 1} to {sign}{end}bp"
            else:
                return f"{sign}{repeat_size}"

    return bin_repeat_size_bp


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
        ~df[f"DiffFromRefRepeats: Allele: Truth"].isna() &
        ~df[f"DiffFromRefRepeats: Allele: {tool}"].isna() & (
                df[f"DiffFromRefRepeats: Allele: {tool}"] != 0) & (
                np.sign(df[f"DiffFromRefRepeats: Allele: Truth"]) != np.sign(df[f"DiffFromRefRepeats: Allele: {tool}"])
        ),
        WRONG_DIRECTION_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        ~df[f"IsRef: Allele: {tool}"].isna() & df[f"IsRef: Allele: {tool}"],
        HET_REF_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        ~df[f"IsHomRef: {tool}"].isna() & df[f"IsHomRef: {tool}"],
        HOM_REF_LABEL,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])


def plot_distribution_by_motif_size(df, figure_title, output_image_path):
    x_column = "MotifSize"
    n_rows = 4
    n_columns = 2
    fig, all_axes = plt.subplots(n_rows, n_columns, figsize=(n_columns*8, n_rows*6), dpi=80, sharex="col", sharey="row")
    for row_i, axes in enumerate(all_axes):
        if row_i == 0:
            hue_column = None
        elif row_i == 1:
            hue_column = "Allele: Concordance: ExpansionHunter vs Truth"
        elif row_i == 2:
            hue_column = "Allele: Concordance: GangSTR vs Truth"
        elif row_i == 3:
            hue_column = "Allele: Concordance: HipSTR vs Truth"
        else:
            raise ValueError(f"Unexpected row_i: {row_i}")

        for column_j, ax in enumerate(axes):
            ax.xaxis.labelpad = ax.yaxis.labelpad = 15
            ax.spines.right.set_visible(False)
            ax.spines.top.set_visible(False)
            ax.tick_params(axis="x")

            ax.set_xlabel("Motif Size (bp)", fontsize=16)
            if hue_column is not None and int(column_j) == 0:
                ax.set_ylabel(hue_column.replace(" vs Truth", "").split(": ")[-1], fontsize=16)

            df_current = df
            ax.set_xlim(left=min(df_current[x_column])-0.5, right=max(df_current[x_column])+0.5)

            title = "Repeats"
            if int(column_j) == 0:
                title = f"Pure " + title
                df_current = df_current[df_current["IsPureRepeat"]]
            elif int(column_j) == 1:
                title += f" With Interruptions"
                df_current = df_current[~df_current["IsPureRepeat"]]
            else:
                raise ValueError(f"Unexpected column_j: {column_j}")

            if row_i != 0:
                title = ""

            p = sns.histplot(
                df_current,
                x=x_column,
                hue=hue_column,
                hue_order=["ExactlyTheSame", "OverlappingCIs", "Discordant"],
                binwidth=1,
                palette=[GREEN_COLOR, "orange", "#FF665533"], # "#5588FF33"
                multiple="stack" if row_i == 0 else "fill",
                stat="proportion",
                discrete=True,
                legend=(column_j == 0) and (row_i > 0),
                ax=axes[column_j])

            p.set_title(title, fontsize=16)

            l = ax.get_legend()
            if l:
                l.set_title("")
                for t, label in zip(l.texts, [
                    "Allele Size Exactly Matches Truth",
                    "Confidence Interval Overlaps Truth",
                    "Confidence Interval Doesn't Overlap Truth",
                ]):
                    t.set_text(label)

    fig.tight_layout()

    print(figure_title)
    suptitle_artist = fig.suptitle(figure_title, fontsize=20, y=1.02)

    plt.savefig(f"{output_image_path}", bbox_extra_artists=(suptitle_artist,), bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_path}")


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
        axes[i].set_xticks(axes[i].get_xticks())
        axes[i].set_xticklabels(
            axes[i].get_xticklabels(),
            rotation=45,
            horizontalalignment="right",
            rotation_mode="anchor")

        if i == 1:
            # add n=.. above each bar
            n_lookup = dict(df.groupby(x_column).count().LocusId)
            for j, (xtick, text) in enumerate(zip(axes[i].get_xticks(), axes[i].get_xticklabels())):
                ax.text(xtick, 1.018, f"{n_lookup[text.get_text()]:,d} alleles", ha="left", va="bottom", color="#777777", rotation=45)

    if tool_name:
        l = axes[0].get_legend()
        l.set_title(f"{tool_name} Call\nvs\nTrue Allele Size\n")
        plt.setp(l.get_title(), multialignment="center")
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


def generate_all_distribution_by_motif_size_plots(df, output_image_dir, plot_counter, max_plots=None, verbose=False):
    for coverage in "40x", "30x", "20x", "10x", "05x", "exome":
        if max_plots and plot_counter >= max_plots:
            print(f"Exiting after generating {plot_counter} plot(s)")
            sys.exit(0)

        df_current = df[
            (df["PositiveOrNegative"] == "positive") &
            (df["coverage"] == coverage) &
            (df["MotifSize"] >= 2) &
            (df["MotifSize"] <= 24) &
            df.IsFoundInReference
        ]

        coverage_label = f"Exome Data" if coverage == "exome" else f"{coverage} Coverage WGS Data"

        output_image_filename = "tool_accuracy_by_motif_size"
        output_image_filename += ".pure_repeats"
        output_image_filename += f".{coverage}"

        output_image_path = f"{output_image_dir}/{output_image_filename}.svg"
        df_current = df_current.sort_values("RepeatSize (bp): Allele: Truth")
        plot_distribution_by_motif_size(df_current, f"Accuracy by Motif Size ({coverage_label})", output_image_path)

        plot_counter += 1

    return plot_counter


def compute_fraction_exactly_correct(df, tool):
    diff_values_considered_exactly_correct = ["0",] # "1", "-1", "2", "-2"]

    exact_match_counts = df[df[f"DiffRepeats: Allele: {tool} - Truth (bin)"].isin(
        diff_values_considered_exactly_correct)].groupby("DiffFromRefRepeats: Allele: Truth (bin)").count()[["LocusId"]]
    total_counts = df.groupby("DiffFromRefRepeats: Allele: Truth (bin)").count()[["LocusId"]]

    total_exact_match = sum(exact_match_counts.LocusId)
    total = sum(total_counts.LocusId)
    overall_fraction_exactly_correct = total_exact_match/total

    result_df = (exact_match_counts/total_counts).fillna(0).rename(columns={"LocusId": "FractionExactlyCorrect"})

    return result_df, overall_fraction_exactly_correct


def compute_tables_for_fraction_exactly_right_plots(df, coverage_values=("40x", "20x", "10x",)):
    tables_by_coverage = []
    for coverage in coverage_values:
        for tool in ("ExpansionHunter", "GangSTR", "HipSTR",):

            df2 = df.copy()
            df2 = df2[df2["coverage"] == coverage]

            df_tool, overall_fraction_exactly_correct = compute_fraction_exactly_correct(df2, tool)
            df_tool.loc[:, "tool"] = f"{tool}: {coverage} coverage" if len(coverage_values) > 1 else tool
            tables_by_coverage.append(df_tool)

            print(f"Processed {tool:20s} --  {100*overall_fraction_exactly_correct:0.1f}% of calls by {tool} "
                  f"@ {coverage} coverage were exactly correct for "
                  f"{len(df2):,d} alleles at {len(set(df2.LocusId)):,d} loci")

    return pd.concat(tables_by_coverage, axis=0)


def generate_all_distribution_by_num_repeats_plots(df, output_image_dir, plot_counter, max_plots=None, verbose=False):

    for motif_size in ("STR", "TR", "TR2"):
        for pure_repeats in (True, False,):
            for coverage in ("40x", "30x", "20x", "10x", "05x", "exome", ):
                for tool_label in ("ExpansionHunter", "GangSTR", "HipSTR"):  #"GangSTR__Q_over_0.8", "ExpansionHunter_Filtered":
                    if max_plots and plot_counter >= max_plots:
                        print(f"Exiting after generating {plot_counter} plot(s)")
                        sys.exit(0)

                    print("-"*100)

                    figure_title_line1 = ""
                    figure_title_line2 = ""
                    output_image_filename = "tool_accuracy_by_true_allele_size"

                    df2 = df.copy()
                    df2 = df2[
                        (df2["PositiveOrNegative"] == "positive") &
                        (df2["coverage"] == coverage) &
                        (df2["DiffFromRefRepeats: Allele: Truth (bin)"] != "0") &
                        df2.IsFoundInReference
                    ]

                    if coverage == "exome":
                        df2 = df2[~df2["Genotype: GangSTR"].isna() | ~df2["Genotype: ExpansionHunter"].isna()]
                        df2 = df2[~df2["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

                    filter_description = []
                    if motif_size == "STR":
                        filter_description.append("2bp to 6bp motifs")
                        output_image_filename += ".2to6bp_motifs"
                        df2 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 6)]
                    elif motif_size == "TR":
                        filter_description.append(f"7bp to 24bp motifs")
                        output_image_filename += ".7to24bp_motifs"
                        df2 = df2[(7 <= df2["MotifSize"]) & (df2["MotifSize"] <= 24)]
                    elif motif_size == "TR2":
                        filter_description.append(f"25bp to 50bp motifs")
                        output_image_filename += ".25to50bp_motifs"
                        df2 = df2[(25 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
                    else:
                        raise ValueError(f"Unexpected motif_size value: {motif_size}")

                    if pure_repeats:
                        #filter_description.append("pure repeats only")
                        output_image_filename += ".pure_repeats"
                        df2 = df2[df2["IsPureRepeat"]]
                    else:
                        filter_description.append("repeats with interruptions")
                        output_image_filename += ".with_interruptions"
                        df2 = df2[~df2["IsPureRepeat"]]

                    output_image_filename += f".{coverage}"

                    tool = tool_label
                    coverage_label = f"Exome Data" if coverage == "exome" else f"{coverage} Coverage Genome Data"
                    if "GangSTR__Q_over_0.8" in tool_label:
                        tool = "GangSTR"
                        figure_title_line1 += f"{tool} Calls filtered to Q > 0.8"
                        output_image_filename += f".{tool}_filtered"
                    elif "ExpansionHunter_Filtered" in tool_label:
                        tool = "ExpansionHunter"
                        figure_title_line1 += f"{tool} Calls filtered"
                        output_image_filename += f".{tool}_filtered"
                    else:
                        figure_title_line1 += f"{tool} Calls for {coverage_label}"
                        output_image_filename += f".{tool}"

                    figure_title_line2 += f"{len(df2):,d} total alleles at {len(set(df2.LocusId)):,d} loci (" + ", ".join(filter_description) + ")"

                    print(figure_title_line1)
                    print(figure_title_line2)

                    skip_condition1 = tool_label == "HipSTR" and motif_size != "STR"
                    skip_condition2 = len(df2) < 10
                    if skip_condition1 or skip_condition2:
                        # HipSTR doesn't support motifs larger than 9bp
                        if skip_condition1:
                            message = "HipSTR doesn't support motif sizes larger than 9bp"
                        elif skip_condition2:
                            message = "Not enough alleles"

                        print(f"Skipping..  {message}")
                        plot_empty_image(
                            figure_title_line1 + "\n\n" + figure_title_line2,
                        )
                        plt.savefig(f"{output_image_dir}/{output_image_filename}.svg")
                        print(f"Saved {output_image_dir}/{output_image_filename}.svg")
                        plt.close()
                        continue

                    print("Keeping", len(df2), "out of ", len(df), "rows")

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
                        palette.append("#880000")
                    if HET_REF_LABEL in hue_values:
                        palette.append("#823d3d")
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
                            print(f"Saved {output_image_dir}/{output_image_filename}{ext}")
                        plt.close()

                    except Exception as e:
                        print(f"ERROR: {e}")
                        #import traceback
                        #traceback.print_exc()

                    plot_counter += 1

    return plot_counter


def generate_fraction_exactly_right_plot(
        df, output_image_dir, plot_counter, discard_hipstr_no_call_loci=False, max_plots=None, verbose=False):
    hue_column = "tool"
    motif_size = "STR"
    pure_repeats = True
    x_column = "DiffFromRefRepeats: Allele: Truth (bin)"
    df = df[
        (df["PositiveOrNegative"] == "positive") &
        (df[x_column] != "0") &
        df.IsFoundInReference
    ]

    if pure_repeats:
        df = df[df["IsPureRepeat"]]
    else:
        df = df[~df["IsPureRepeat"]]

    if motif_size == "STR":
        df = df[(2 <= df["MotifSize"]) & (df["MotifSize"] <= 6)]
    elif motif_size == "TR":
        df = df[(7 <= df["MotifSize"]) & (df["MotifSize"] <= 24)]
    elif motif_size == "TR2":
        df = df[(25 <= df["MotifSize"]) & (df["MotifSize"] <= 50)]
    else:
        raise ValueError(f"Unexpected motif_size value: {motif_size}")

    if discard_hipstr_no_call_loci:
        count_before = len(set(df[df["coverage"] == "40x"].LocusId))
        df = df[df[f"DiffRepeats: Allele: HipSTR - Truth (bin)"] != NO_CALL_LABEL]
        count_discarded = count_before - len(set(df[df["coverage"] == "40x"].LocusId))
        print(f"Discarded {count_discarded:,d} ({100.0*count_discarded/count_before:0.1f}%) of loci due to HipSTR no call")

    df_fraction = compute_tables_for_fraction_exactly_right_plots(df, coverage_values=("40x", "30x", "20x", "10x", "05x"))

    df_fraction = df_fraction.reset_index()
    df_fraction = df_fraction.sort_values(x_column, key=lambda c: c.str.split(" ").str[0].astype(int))

    for coverage in "40x", "all":
        df2 = df_fraction.copy()
        if coverage == "all":
            # exclude HipSTR to reduce clutter, and also 30x and 5x coverage
            df2 = df2[
                ~df2.tool.str.contains("HipSTR") &
                ~df2.tool.str.contains("30x") &
                ~df2.tool.str.contains("05x")
            ]
        else:
            df2 = df2[df2.tool.str.contains(coverage)]
            df2.loc[:, "tool"] = df2["tool"].apply(lambda s: s.split(":")[0])

        filter_description = []
        filename_suffix = ""
        if motif_size == "STR":
            filter_description.append("2bp to 6bp motifs")
            #figure_title += "2bp to 6bp motifs"
            filename_suffix += ".2to6bp_motifs"
        elif motif_size == "TR":
            filter_description.append(f"7bp to 24bp motifs")
            #figure_title += f"7bp to 24bp motifs"
            filename_suffix += ".7to24bp_motifs"
        elif motif_size == "TR2":
            filter_description.append(f"25bp to 50bp motifs")
            #figure_title += f"7bp to 24bp motifs"
            filename_suffix += ".25to50bp_motifs"
        else:
            raise ValueError(f"Unexpected motif_size value: {motif_size}")

        if pure_repeats:
            #filter_description.append("pure repeats only")
            #figure_title += ", pure repeats only"
            filename_suffix += ".pure_repeats"
        else:
            filter_description.append("repeats with interruptions")
            #figure_title += "only repeats with interruptions"
            filename_suffix += ".has_interruptions"

        if coverage == "all":
            filename_suffix += ".compare_different_coverages"
        else:
            filename_suffix += f".{coverage}_coverage"

        if discard_hipstr_no_call_loci:
            filter_description.append("excluding HipSTR no-call loci")
            filename_suffix += f".excluding_hipstr_no_call_loci"

        figure_title = f"Accuracy of " + ", ".join(sorted(set([t.split(":")[0] for t in set(df2.tool)])))
        figure_title += "\n\n"
        figure_title += f"at {len(set(df.LocusId)):,d} {motif_size.strip('2')} loci (" + ", ".join(filter_description) + ")"

        print(figure_title.replace("\n", " "))

        hue_values = set(df2[hue_column])
        num_gangstr_colors = len([h for h in hue_values if h.lower().startswith("gangstr")])
        num_expansion_hunter_colors = len([h for h in hue_values if h.lower().startswith("expansionhunter")])
        num_hipstr_colors = len([h for h in hue_values if h.lower().startswith("hipstr")])

        # plot figure
        fig, ax = plt.subplots(1, 1, figsize=(12, 10), dpi=80)

        def hue_order(h):
            tokens = h.split(": ")
            try:
                return tokens[0], -1*int(tokens[1][0:2])
            except Exception as e:
                if verbose:
                    print(f"Unable to parse hue value: {h}: {e}")
                return h

        sns.pointplot(
            data=df2,
            x=x_column,
            y="FractionExactlyCorrect",
            hue=hue_column,
            scale=0.8,
            hue_order=sorted(hue_values, key=hue_order),
            palette=(
                    list(sns.color_palette("Purples_r", n_colors=num_expansion_hunter_colors) if num_expansion_hunter_colors > 1 else ["#6A51A3"]) +
                    list(sns.color_palette("blend:#FF3355,#FFCCCC", n_colors=num_gangstr_colors)) +
                    list(sns.color_palette("Oranges_r", n_colors=num_hipstr_colors))
            ),
            ax=ax,
        )

        ax.set_xlabel("True Allele Size Minus Number of Repeats in Reference Genome", fontsize=16)
        ax.set_ylabel("Fraction of Calls That Exactly Match True Allele Size", fontsize=16)

        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.set_ylim((0, 1))
        ax.yaxis.set_ticks(np.arange(0, 1.05, 0.1))
        ax.grid(axis='y', color='#ECECEC')

        ax.set_xticks(ax.get_xticks())
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment="right", rotation_mode='anchor')
        ax.get_legend().set_title(f"")
        ax.legend(frameon=True)
        if len(hue_values) > 4:
            ax.get_legend().set_bbox_to_anchor((0.15, 0.25))
        else:
            #ax.legend(loc="upper right", frameon=True)
            ax.get_legend().set_bbox_to_anchor((0.15, 0.15))

        #fig.tight_layout()
        suptitle_artist = fig.suptitle(figure_title, fontsize=17, y=1.01)

        output_image_filename = "tool_accuracy_by_true_allele_size_exactly_matching_calls"

        for ext in ".svg", ".png":
            plt.savefig(f"{output_image_dir}/{output_image_filename}{filename_suffix}{ext}", bbox_extra_artists=(suptitle_artist,), bbox_inches="tight")
            print(f"Saved {output_image_dir}/{output_image_filename}{filename_suffix}{ext}")

        plt.close()
        plot_counter += 1

        if max_plots and plot_counter >= max_plots:
            print(f"Exiting after generating {plot_counter} plot(s)")
            sys.exit(0)

    return plot_counter


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="If specified, only generate this many plots. Useful for testing")
    p.add_argument("--skip-plot1", action="store_true", help="Don't generate line plots with fraction exactly right")
    p.add_argument("--skip-plot2", action="store_true", help="Don't generate bar plots by motif size")
    p.add_argument("--skip-plot3", action="store_true", help="Don't generate bar plots by true allele size")
    p.add_argument("--output-dir", default="../tool_comparison/figures")
    p.add_argument("--verbose", action="store_true", help="Print additional info")
    args = p.parse_args()

    input_table_path = "../tool_comparison/combined.results.alleles.tsv"

    print(f"Loading {input_table_path}")
    df = pd.read_table(input_table_path)

    #df = df[df.Motif.isin({"AAG", "AGA", "GAA", "CTT", "TCT", "TTC"})]
    #df = df[df.MotifSize == 3]

    if args.verbose:
        print("Num loci:")
        print(df.groupby(["PositiveOrNegative", "coverage", "IsPureRepeat"]).count().LocusId/2)

    print("Computing additional columns...")
    df.loc[:, "RepeatSize (bp): Allele: Truth (bin)"] = df["RepeatSize (bp): Allele: Truth"].apply(
        bin_repeat_size_bp_wrapper(400, 24))

    df.loc[:, "DiffFromRefRepeats: Allele: Truth (bin)"] = df.apply(bin_num_repeats_wrapper(bin_size=2), axis=1)
    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    for tool in ("ExpansionHunter", "GangSTR", "HipSTR",):
        define_hue_column(df, tool)

    print("Generating plots...")
    plot_counter = 0
    if not args.skip_plot1:
        plot_counter = generate_fraction_exactly_right_plot(
            df, args.output_dir, plot_counter, discard_hipstr_no_call_loci=False, max_plots=args.n, verbose=args.verbose)
        plot_counter = generate_fraction_exactly_right_plot(
            df, args.output_dir, plot_counter, discard_hipstr_no_call_loci=True, max_plots=args.n, verbose=args.verbose)

    if not args.skip_plot2:
        plot_counter = generate_all_distribution_by_motif_size_plots(
            df, args.output_dir, plot_counter, max_plots=args.n, verbose=args.verbose)
    if not args.skip_plot3:
        plot_counter = generate_all_distribution_by_num_repeats_plots(
            df, args.output_dir, plot_counter, max_plots=args.n, verbose=args.verbose)

    print(f"Done generating all {plot_counter} plots")


if __name__ == "__main__":
    main()