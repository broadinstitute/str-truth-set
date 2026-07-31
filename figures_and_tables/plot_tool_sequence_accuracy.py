"""Plot how close a tool's reported allele *sequences* are to the assembly truth allele sequences.

This is the sequence-accuracy counterpart of plot_tool_accuracy_by_allele_size.py. That script scores a tool by how
close its repeat count is to the truth; this one scores the tools that report an actual allele sequence (TRGT,
ATaRVa, HipSTR) by the edit distance columns that add_sequence_accuracy_columns.py adds:

    SequenceEditDistance: Allele: {tool}             Levenshtein distance in bp
    SequenceEditDistanceNormalized: Allele: {tool}   the same distance / max(1, len(truth allele sequence))

Each figure has two panels sharing the accuracy figure's x axis (true allele size minus the reference repeat count):

    left    a box plot of the edit distance per x bin, so the median and spread are readable directly in bp
    right   a 100%-stacked histogram of edit-distance categories per x bin, so "how much of this bin is
            sequence-exact" reads at a glance

The left panel deliberately differs from the accuracy figure's left panel, which is a stacked count: a count of
alleles says nothing about how far off the sequences are, which is the whole point here.
"""

import argparse
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns

from plot_tool_accuracy_by_allele_size import (
    GREEN_COLOR, NO_CALL_LABEL, SEQUENCING_DATA_TYPE_LABELS, TITLE_TOOL_LABELS,
    bin_num_repeats_wrapper, plot_empty_image)

sns.set_context(font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
})

# Tools scored by this benchmark. Must match run_tools/run_genotyping_tools.py's SEQUENCE_ACCURACY_TOOLS --
# TRGTv3 is scoped out of v1 there, so it's excluded here too even though its VCF format is TRGTv5-compatible.
SEQUENCE_ACCURACY_TOOLS = ["TRGTv5", "ATaRVa", "HipSTR"]

# (upper bound inclusive, label) for the stacked category panel, coarsest last. The bp cuts are absolute; the
# normalized cuts are fractions of the true allele length.
EDIT_DISTANCE_BINS = [
    (0, "0 (exact)"),
    (2, "1 to 2 bp"),
    (5, "3 to 5 bp"),
    (10, "6 to 10 bp"),
    (float("inf"), "more than 10 bp"),
]

NORMALIZED_EDIT_DISTANCE_BINS = [
    (0, "0 (exact)"),
    (0.01, "up to 1%"),
    (0.05, "up to 5%"),
    (0.10, "up to 10%"),
    (float("inf"), "more than 10%"),
]

# green -> red, in the same order as the *_EDIT_DISTANCE_BINS labels above
CATEGORY_COLORS = [GREEN_COLOR, "#A3CC53", "#F2C744", "#EE8B3A", "#C9342A"]
NO_CALL_COLOR = "#c5c5c5"

# --metric token -> the column it reads, the category cuts it uses, and its y axis label. The token is also the
# ".{metric}" component of the output filename, which docs/tool_comparison_viewer.html builds its image URLs from.
METRICS = {
    "edit_distance_bp": {
        "column_prefix": "SequenceEditDistance",
        "bins": EDIT_DISTANCE_BINS,
        "axis_label": "Edit Distance vs. True Allele Sequence (bp)",
    },
    "edit_distance_normalized": {
        "column_prefix": "SequenceEditDistanceNormalized",
        "bins": NORMALIZED_EDIT_DISTANCE_BINS,
        "axis_label": "Edit Distance / True Allele Length",
    },
}

X_COLUMN = "DiffFromRefRepeats: Allele: Truth (bin)"


def bin_edit_distance(distance, bins):
    """Convert one edit distance to its category label.

    Args:
        distance (float): the edit distance, or NaN when the allele couldn't be scored.
        bins (list): (upper bound inclusive, label) tuples, coarsest last.

    Returns:
        str: the category label, or NO_CALL_LABEL when the distance is NaN.
    """
    if pd.isna(distance):
        return NO_CALL_LABEL
    for upper_bound, label in bins:
        if distance <= upper_bound:
            return label
    raise ValueError(f"Unexpected edit distance: {distance}")


def plot_sequence_accuracy(df, args, distance_column, category_column, category_order, palette, axis_label,
                           tool_name=None, figure_title=None):
    """Draw the two-panel figure: a box plot of the distances on the left, stacked categories on the right."""
    fig, axes = plt.subplots(1, 2, figsize=(args.width, args.height))

    if figure_title:
        fig.suptitle(figure_title, fontsize=17)

    # the x bins are strings, so give both panels the same explicit order (numeric order of the underlying diff)
    x_order = list(dict.fromkeys(df[X_COLUMN]))

    for i, ax in enumerate(axes):
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.set_xlabel("True Allele Size Minus Number of Repeats in Reference Genome", fontsize=14)
        ax.tick_params(axis="x")

    axes[0].spines.right.set_visible(False)
    axes[0].spines.top.set_visible(False)
    axes[0].set_ylabel(axis_label, fontsize=14)
    sns.boxplot(
        data=df[df[distance_column].notna()],
        x=X_COLUMN,
        y=distance_column,
        order=x_order,
        color=GREEN_COLOR,
        showfliers=False,
        ax=axes[0])

    axes[1].set_ylabel("Fraction of Alleles", fontsize=14)
    sns.histplot(
        df,
        x=X_COLUMN,
        hue=category_column,
        hue_order=category_order,
        binwidth=1,
        palette=palette,
        multiple="fill",
        stat="proportion",
        discrete=True,
        legend=True,
        ax=axes[1])

    fig.tight_layout()

    n_lookup = dict(df.groupby(X_COLUMN).count().LocusId)
    for ax in axes:
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment="right", rotation_mode="anchor",
                           fontsize=11)
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=12)
        # add n=.. above each bar, as in plot_tool_accuracy_by_allele_size
        for xtick, text in zip(ax.get_xticks(), ax.get_xticklabels()):
            n_for_bin = n_lookup.get(text.get_text())
            if n_for_bin:
                ax.text(xtick, 1.018, f"{n_for_bin:,d}", ha="left", va="bottom", color="#777777", fontsize=12,
                        rotation=45, transform=ax.get_xaxis_transform())
        ax.text(0, 1.15, "Alleles Per Bin", ha="left", va="bottom", color="#777777", fontsize=12,
                transform=ax.transAxes)

    if tool_name:
        legend = axes[1].get_legend()
        legend.set_title(f"{tool_name} Allele Sequence\nvs\nTrue Allele Sequence\n")
        legend.get_title().set_horizontalalignment("center")
        legend.set_frame_on(False)

    fig.tight_layout()


def generate_all_plots(df, args):
    """Generate one figure per (metric x motif size bin x genotype subset) and save it to args.output_dir.

    Returns:
        int: the number of figures written.
    """
    plot_counter = 0
    for metric_name in [args.metric] if args.metric else list(METRICS):
        metric = METRICS[metric_name]
        distance_column = f"{metric['column_prefix']}: Allele: {args.tool}"
        if distance_column not in df.columns:
            print(f"WARNING: {distance_column} not found in the input table. Skipping {metric_name}...")
            continue

        category_column = f"{distance_column} (bin)"
        df.loc[:, category_column] = [bin_edit_distance(d, metric["bins"]) for d in df[distance_column]]

        for motif_size in (["all_motifs"] if args.min_motif_size is None and args.max_motif_size is None
                           else [f"{args.min_motif_size}-{args.max_motif_size}bp"]):
            if motif_size == "all_motifs":
                df_motif = df[(2 <= df["MotifSize"]) & (df["MotifSize"] <= 50)]
                motif_filename_token = ".all_motifs"
                motif_description = "all motif sizes"
            else:
                df_motif = df
                if args.min_motif_size:
                    df_motif = df_motif[args.min_motif_size <= df_motif["MotifSize"]]
                if args.max_motif_size:
                    df_motif = df_motif[df_motif["MotifSize"] <= args.max_motif_size]
                motif_filename_token = f".{args.min_motif_size}to{args.max_motif_size}bp_motifs"
                motif_description = f"{args.min_motif_size}bp to {args.max_motif_size}bp motifs"

            for genotype_subset in ["all", "HET", "HOM"] if not args.genotype else [args.genotype]:
                if genotype_subset == "all":
                    df_plot = df_motif
                    genotype_filename_token = ".all_genotypes"
                else:
                    df_plot = df_motif[df_motif["SummaryString"].str.contains(f":{genotype_subset}")]
                    genotype_filename_token = f".{genotype_subset}"

                output_image_filename = ("tool_sequence_accuracy_by_true_allele_size"
                                         f".{metric_name}{motif_filename_token}{genotype_filename_token}"
                                         f".{args.coverage}.{args.tool}")

                filter_description = [motif_description]
                if genotype_subset != "all":
                    filter_description.append(genotype_subset)

                data_type_label = SEQUENCING_DATA_TYPE_LABELS.get(args.sequencing_data_type, "genome")
                coverage_label = ("exome data" if args.coverage == "exome"
                                  else f"{args.coverage} {data_type_label} data")

                # HipSTR doesn't support motifs larger than 9bp, so bins entirely above that are always empty.
                # A None lower bound means "unstratified" (the all_motifs bin), not "no lower limit" -- it must
                # not be treated as satisfying "> 9", or the all_motifs bin (which does include HipSTR's real
                # 2-9bp calls) gets wrongly blanked too.
                skip_message = None
                if args.tool == "HipSTR" and (
                        args.min_motif_size is not None and args.min_motif_size > 9
                        and (args.max_motif_size is None or args.max_motif_size > 9)):
                    skip_message = "HipSTR doesn't support motif sizes larger than 9bp"
                elif len(df_plot) < 10:
                    skip_message = "Not enough alleles to create plot"

                n_locus_ids = len(set(df_plot.LocusId))
                figure_title_line2 = f"{n_locus_ids:,d} loci (" + ", ".join(filter_description) + ")"
                if skip_message:
                    figure_title_line1 = f"{TITLE_TOOL_LABELS.get(args.tool, args.tool)} {coverage_label}"
                    print(f"Skipping {output_image_filename}: {skip_message}")
                    plot_empty_image(figure_title_line1 + "\n\n" + figure_title_line2, skip_message,
                                     args.width, args.height)
                else:
                    num_alleles_exactly_right = sum(df_plot[category_column] == metric["bins"][0][1])
                    figure_title_line1 = (
                        f"{TITLE_TOOL_LABELS.get(args.tool, args.tool)} got {num_alleles_exactly_right:,d} out of "
                        f"{len(df_plot):,d} allele sequences ({100*num_alleles_exactly_right/len(df_plot):0.1f}%) "
                        f"exactly right in {coverage_label}")
                    print(figure_title_line1)
                    print(f"Plotting {len(df_plot):,d} out of {len(df):,d} rows")

                    category_order = [NO_CALL_LABEL] + [label for _, label in metric["bins"]]
                    present = set(df_plot[category_column])
                    category_order = [label for label in category_order if label in present]
                    palette = [NO_CALL_COLOR] + CATEGORY_COLORS
                    palette = [color for label, color in zip(
                        [NO_CALL_LABEL] + [l for _, l in metric["bins"]], palette) if label in present]

                    plot_sequence_accuracy(
                        df_plot, args,
                        distance_column=distance_column,
                        category_column=category_column,
                        category_order=category_order,
                        palette=palette,
                        axis_label=metric["axis_label"],
                        tool_name=args.tool,
                        figure_title=figure_title_line1 + "\n\n" + figure_title_line2 if args.show_title else None)

                output_path = os.path.join(args.output_dir, f"{output_image_filename}.{args.image_type}")
                plt.savefig(output_path, dpi=300)
                plt.close()
                print(f"Saved {output_path}")
                plot_counter += 1

    return plot_counter


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--width", default=20, type=float)
    p.add_argument("--height", default=9, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"], help="Image type to generate")
    p.add_argument("--show-title", action="store_true", help="Show plot title")
    p.add_argument("--output-dir", default=".")

    g = p.add_argument_group("Filters")
    g.add_argument("--tool", choices=SEQUENCE_ACCURACY_TOOLS, required=True, help="Which tool to plot")
    g.add_argument("--metric", choices=sorted(METRICS), help="Plot only this metric. Both are plotted by default.")
    g.add_argument("--coverage", required=True, help="Coverage label of the input table (example: \"30x\")")
    p.add_argument("--sequencing-data-type", choices=sorted(SEQUENCING_DATA_TYPE_LABELS),
                   help="Sequencing data type, used to label the plot title (example: \"pacbio\" -> \"PacBio HiFi\").")
    g.add_argument("--min-motif-size", type=int, help="Min motif size. Without --min/--max-motif-size the script "
                   "generates the unstratified all_motifs plot (so it needs no --all-motifs-only flag).")
    g.add_argument("--max-motif-size", type=int, help="Max motif size")
    g.add_argument("--genotype", choices=["all", "HET", "HOM"], help="Plot only this genotype subset")

    p.add_argument("combined_tool_results_tsv", help="Path of the ....alleles.tsv.gz table")
    args = p.parse_args()

    print("Command line args:")
    for arg in vars(args):
        print(f"{arg}: {getattr(args, arg)}")

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)

    before = len(df)
    df = df[~df["DiffFromRefRepeats: Allele: Truth"].isna()]
    print(f"Filtered out {before - len(df):,d} of {before:,d} rows with no true allele size")

    if "Coverage" not in df.columns:
        df["Coverage"] = args.coverage
    df = df[df["Coverage"] == args.coverage]
    print(f"Kept {len(df):,d} rows at coverage {args.coverage}")

    if args.coverage == "exome" and "GeneRegionFromGencode_V42" in df.columns:
        df = df[~df["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

    if len(df) == 0:
        raise ValueError(f"No rows left in {args.combined_tool_results_tsv} after filtering to "
                         f"coverage {args.coverage}")

    print("Computing additional columns...")
    df.loc[:, X_COLUMN] = df.apply(bin_num_repeats_wrapper(bin_size=2), axis=1)
    df = df.sort_values("DiffFromRefRepeats: Allele: Truth")

    print(f"Generated {generate_all_plots(df, args):,d} plots")


if __name__ == "__main__":
    main()
