import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
    "legend.fontsize": 14,
    "legend.loc": "upper right",
})


GREEN_COLOR = "#50AA44"

NO_CALL_LABEL = "No Call"


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


def compute_fraction_exactly_right(df, tool):
    subset_columns = ["DiffFromRefRepeats: Allele: Truth (bin)", "LocusId"]

    total_counts = df[subset_columns].groupby("DiffFromRefRepeats: Allele: Truth (bin)").count()

    df_exactly_right = df[df[f"DiffRepeats: Allele: {tool} - Truth"] == 0]
    exact_match_counts = df_exactly_right[subset_columns].groupby("DiffFromRefRepeats: Allele: Truth (bin)").count()

    total_exact_match = sum(exact_match_counts.LocusId)
    total = sum(total_counts.LocusId)
    overall_fraction_exactly_right = total_exact_match/total

    result_df = (exact_match_counts/total_counts).fillna(0).rename(columns={"LocusId": "FractionExactlyRight"})

    return result_df, overall_fraction_exactly_right


def compute_tables_for_fraction_exactly_right_plots(df, coverage_values=("40x", "20x", "10x",)):
    tables_by_coverage = []
    for coverage in coverage_values:
        for tool in ("ExpansionHunter", "GangSTR", "HipSTR"):

            df_current = df.copy()
            df_current = df_current[df_current["Coverage"] == coverage]

            df_tool, overall_fraction_exactly_right = compute_fraction_exactly_right(df_current, tool)
            df_tool.loc[:, "tool"] = f"{tool}: {coverage} coverage" if len(coverage_values) > 1 else tool
            tables_by_coverage.append(df_tool)

            print(f"Processed {tool:20s} --  {100*overall_fraction_exactly_right:0.1f}% of calls by {tool} "
                  f"@ {coverage} coverage were exactly right for "
                  f"{len(df_current):,d} alleles at {len(set(df_current.LocusId)):,d} loci")

    return pd.concat(tables_by_coverage, axis=0)


def generate_fraction_exactly_right_plot(df, output_dir, exclude_hipstr_no_call_loci=False, verbose=False, show_title=True):

    df_fraction = compute_tables_for_fraction_exactly_right_plots(df, coverage_values=("40x","20x", "10x"))
    df_fraction = df_fraction.reset_index()
    df_fraction = df_fraction.sort_values("DiffFromRefRepeats: Allele: Truth (bin)", key=lambda c: c.str.split(" ").str[0].astype(int))

    for coverage in "40x", "all":
        df_current = df_fraction.copy()
        if coverage == "all":
            # exclude HipSTR to make plot clearer
            df_current = df_current[~df_current.tool.str.contains("HipSTR")]
        else:
            df_current = df_current[df_current.tool.str.contains(coverage)]
            df_current.loc[:, "tool"] = df_current["tool"].apply(lambda s: s.split(":")[0])

        filter_description = ["2bp to 6bp motifs"]
        filename_suffix = ".2to6bp_motifs.pure_repeats"
        if coverage == "all":
            filename_suffix += ".compare_different_coverages"
        else:
            filename_suffix += f".{coverage}_coverage"

        if exclude_hipstr_no_call_loci:
            filter_description.append("excluding HipSTR no-call loci")
            filename_suffix += f".excluding_hipstr_no_call_loci"

        figure_title = f"Accuracy of " + ", ".join(sorted(set([t.split(":")[0] for t in set(df_current.tool)])))
        figure_title += "\n\n"
        figure_title += f"at {len(set(df.LocusId)):,d} STR loci (" + ", ".join(filter_description) + ")"

        print("Plotting", figure_title.replace("\n", " "))

        hue_values = set(df_current["tool"])
        num_gangstr_colors = len([h for h in hue_values if h.lower().startswith("gangstr")])
        num_expansion_hunter_colors = len([h for h in hue_values if h.lower().startswith("expansionhunter")])
        num_hipstr_colors = len([h for h in hue_values if h.lower().startswith("hipstr")])

        # plot figure
        fig, ax = plt.subplots(1, 1, figsize=(12, 10))

        def hue_order(h):
            tokens = h.split(": ")
            try:
                return tokens[0], -1*int(tokens[1][0:2])
            except Exception as e:
                if verbose:
                    print(f"Unable to parse hue value: {h}: {e}")
                return h

        sns.pointplot(
            data=df_current,
            x="DiffFromRefRepeats: Allele: Truth (bin)",
            y="FractionExactlyRight",
            hue="tool",
            scale=0.8,
            hue_order=sorted(hue_values, key=hue_order),
            palette=(
                list(sns.color_palette("Purples_r", n_colors=num_expansion_hunter_colors)
                     if num_expansion_hunter_colors > 1 else ["#6A51A3"]) +
                list(sns.color_palette("blend:#FF3355,#FFCCCC", n_colors=num_gangstr_colors)) +
                list(sns.color_palette("Oranges_r", n_colors=num_hipstr_colors))
            ),
            ax=ax,
        )

        ax.grid(axis='y', color='#ECECEC')

        ax.set_xlabel("True Allele Size Minus Number of Repeats in Reference Genome", fontsize=16)
        ax.set_ylabel("Fraction of Calls That Exactly Match True Allele Size", fontsize=16)

        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment="right", rotation_mode='anchor', fontsize=14)
        y_ticks = np.arange(0, 1.05, 0.1)
        ax.set_ylim((0, 1))
        ax.set_yticks(y_ticks)
        ax.set_yticklabels([f"{t:0.1f}" for t in y_ticks], fontsize=15)

        ax.get_legend().set_title(f"")
        ax.get_legend().set_frame_on(False)
        if coverage == "all":
            ax.get_legend().set_bbox_to_anchor((0.15, 0.35))

        if show_title:
            suptitle_artist = fig.suptitle(figure_title, fontsize=17, y=1.01)
            extra_artists = [suptitle_artist]
        else:
            extra_artists = []

        output_image_filename = "tool_accuracy_by_true_allele_size_exactly_matching_calls"

        plt.savefig(os.path.join(output_dir, f"{output_image_filename}{filename_suffix}.svg"),
                    bbox_extra_artists=extra_artists, bbox_inches="tight")
        print(f"Saved {output_dir}/{output_image_filename}{filename_suffix}.svg")

        plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--exclude-hipstr-no-call-loci", action="store_true", help="In the plot, exclude HipSTR no-call loci")
    p.add_argument("--output-dir", default=".")
    p.add_argument("--verbose", action="store_true", help="Print additional info")
    p.add_argument("combined_tool_results_tsv", nargs="?", default="../tool_comparison/combined.results.alleles.tsv.gz")
    args = p.parse_args()

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)
    print(f"Loaded {len(df)} from {args.combined_tool_results_tsv}")

    print("Filtering table...")
    df = df[(df["PositiveOrNegative"] == "positive") & df["IsFoundInReference"] & df["IsPureRepeat"] &
            (2 <= df["MotifSize"]) & (df["MotifSize"] <= 6)]

    if args.verbose:
        print("Num loci:")
        print(df.groupby(["PositiveOrNegative", "Coverage", "IsPureRepeat"]).count().LocusId/2)

    if args.exclude_hipstr_no_call_loci:
        count_before = len(set(df[df["Coverage"] == "40x"].LocusId))
        df = df[~df[f"DiffRepeats: Allele: HipSTR - Truth"].isna()]
        count_discarded = count_before - len(set(df[df["Coverage"] == "40x"].LocusId))
        print(f"Discarded {count_discarded:,d} out of {count_before:,d} ({100.0*count_discarded/count_before:0.1f}%) "
              f"of loci due to HipSTR no call")

    print("Computing bin column: 'DiffFromRefRepeats: Allele: Truth (bin)' ")
    df.loc[:, "DiffFromRefRepeats: Allele: Truth (bin)"] = df.apply(bin_num_repeats_wrapper(bin_size=2), axis=1)

    print("Sorting by 'DiffFromRefRepeats: Allele: Truth' column")
    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    generate_fraction_exactly_right_plot(
        df, args.output_dir, exclude_hipstr_no_call_loci=args.exclude_hipstr_no_call_loci, verbose=args.verbose)


if __name__ == "__main__":
    main()