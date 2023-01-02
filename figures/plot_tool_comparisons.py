import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

sns.set_context("notebook", font_scale=1.1, rc={
    "font.family": "sans-serif",
})

GREEN_COLOR = "#50AA44"


def bin_num_repeats_wrapper(threshold=20, bin_size=1):
    """Create a function that converts the integer number of STR repeats to a histogram bin"""

    threshold = threshold - threshold % bin_size + bin_size + 1
    def bin_num_repeats(num_repeats):
        sign = "-" if num_repeats < 0 else "+"
        num_repeats = abs(num_repeats)
        if num_repeats == 0:
            return "0"
        elif num_repeats >= threshold:
            return f"{sign}{int(threshold)}" + " or more"
        else:
            if bin_size > 1:
                start = int((num_repeats - 1)/bin_size) * bin_size
                end = start + bin_size
                return f"{sign}{int(start + 1)} to {sign}{int(end)}"
            else:
                return f"{sign}{int(num_repeats)}"

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
            return "No Call"

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

    wrong_direction_label = "Wrong Direction"
    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        ~df[f"DiffFromRefRepeats: Allele: Truth"].isna() &
        ~df[f"DiffFromRefRepeats: Allele: {tool}"].isna() & (
                df[f"DiffFromRefRepeats: Allele: {tool}"] != 0) & (
                np.sign(df[f"DiffFromRefRepeats: Allele: Truth"]) != np.sign(df[f"DiffFromRefRepeats: Allele: {tool}"])
        ),
        wrong_direction_label,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    ref_allele_label = "Called Het Ref"
    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        ~df[f"IsRef: Allele: {tool}"].isna() & df[f"IsRef: Allele: {tool}"],
        ref_allele_label,
        df[f"DiffRepeats: Allele: {tool} - Truth (bin)"])

    ref_variant_label = "Called Hom Ref"
    df.loc[:, f"DiffRepeats: Allele: {tool} - Truth (bin)"] = np.where(
        ~df[f"IsHomRef: {tool}"].isna() & df[f"IsHomRef: {tool}"],
        ref_variant_label,
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

    print(f"Saved {output_image_path}")
    plt.savefig(f"{output_image_path}", bbox_extra_artists=(suptitle_artist,), bbox_inches="tight")
    plt.close()


def plot_distribution_by_num_repeats(
        df,
        x_column,
        hue_column,
        hue_order=None,
        tool_name=None,
        palette=None,
        figure_title=None,
):

    n_plots = 2
    fig, axes = plt.subplots(1, n_plots, figsize=(n_plots*10, 9), dpi=100)

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
    if value == "Filtered":
        return -2500
    elif value == "No Call":
        return -2000
    elif value == "Called Hom Ref":
        return -1500
    elif value == "Called Het Ref":
        return -1000
    elif value == "Wrong Direction":
        return -500
    else:
        return int(float(value.split(" ")[0]))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="If specified, only generate this many plots. Useful for testing")
    args = p.parse_args()

    input_table_path = "../tool_comparison/combined.results.alleles.tsv"
    output_image_dir = "../tool_comparison/figures"
    print(f"Loading {input_table_path}")
    df = pd.read_table(input_table_path)

    #print("Num loci:")
    #print(df.groupby(["PositiveOrNegative", "coverage", "IsPureRepeat"]).count()["LocusId"]/2)

    df.loc[:, "RepeatSize (bp): Allele: Truth (bin)"] = df["RepeatSize (bp): Allele: Truth"].apply(
        bin_repeat_size_bp_wrapper(400, 24))

    plot_counter = 0
    for coverage in "40x", "30x", "20x", "10x", "05x", "exome":
        if args.n and plot_counter >= args.n:
            print(f"Exiting after generating {args.n} plot(s)")
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

    df.loc[:, "DiffFromRefRepeats: Allele: Truth (bin)"] = df["DiffFromRefRepeats: Allele: Truth"].apply(
        bin_num_repeats_wrapper(19, bin_size=2))
    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    for motif_size in ("STR", "TR"):
        for pure_repeats in (True, False,):
            for coverage in ("40x", "30x", "20x", "10x", "05x", "exome", ):
                for tool_label in ("ExpansionHunter", "GangSTR", "HipSTR"):  #"GangSTR__Q_over_0.8", "ExpansionHunter_Filtered":
                    if tool_label == "HipSTR" and motif_size == "TR":
                        # HipSTR doesn't support motifs larger than 9bp
                        continue

                    if args.n and plot_counter >= args.n:
                        print(f"Exiting after generating {args.n} plot(s)")
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
                    else:
                        raise ValueError(f"Unexpected motif_size value: {motif_size}")

                    if pure_repeats:
                        filter_description.append("pure repeats only")
                        output_image_filename += ".pure_repeats"
                        df2 = df2[df2["IsPureRepeat"]]
                    else:
                        filter_description.append("only repeats with interruptions")
                        output_image_filename += ".with_interruptions"
                        df2 = df2[~df2["IsPureRepeat"]]

                    output_image_filename += f".{coverage}"

                    tool = tool_label
                    coverage_label = f"Exome Data" if coverage == "exome" else f"{coverage} Coverage WGS Data"
                    if "GangSTR__Q_over_0.8" in tool_label:
                        tool = "GangSTR"
                        figure_title_line1 += f"{tool} Calls filtered to Q > 0.8"
                        output_image_filename += f".{tool}_filtered"
                    elif "ExpansionHunter_Filtered" in tool_label:
                        tool = "ExpansionHunter"
                        figure_title_line1 += f"{tool} Calls filtered"
                        output_image_filename += f".{tool}_filtered"
                    else:
                        figure_title_line1 += f"{tool} Calls in {coverage_label}"
                        output_image_filename += f".{tool}"

                    figure_title_line2 += f"Showing {len(df2):,d} total alleles at {len(set(df2.LocusId)):,d} loci (" + ", ".join(filter_description) + ")"

                    print(figure_title_line1)
                    print(figure_title_line2)
                    if len(df2) < 1000:
                        print(f"Skipping..  only {len(df)} total alleles")
                        continue

                    print("Keeping", len(df2), "out of ", len(df), "rows")

                    print(tool, sum(df2[f"DiffFromRefRepeats: Allele: {tool}"].isna()), "out of", len(df2), f"{tool} - Ref' values are NaN")
                    print(tool, sum(df2[f"DiffRepeats: Allele: {tool} - Truth"].isna()), "out of", len(df2), f"{tool} - Truth' values are NaN")

                    define_hue_column(df2, tool)

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

                    palette = ["#c9c9c9", "#880000", "#823d3d", "#99FFEF"]
                    if "Filtered" in hue_values:
                        palette = ["#f5f5f5"] + palette

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

                        print(f"Saved {output_image_dir}/{output_image_filename}.svg")
                        plt.savefig(f"{output_image_dir}/{output_image_filename}.svg")
                        plt.close()

                    except Exception as e:
                        print(f"ERROR: {e}")
                        #import traceback
                        #traceback.print_exc()

                    plot_counter += 1


if __name__ == "__main__":
    main()