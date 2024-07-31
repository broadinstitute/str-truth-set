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

sns.set_palette(["purple", "blue", "red", "orange"])


TITLE_FONT_SIZE = 13


def plot(df, width, height, figure_title=None, only_expansion_hunter=False):
    total_allele_count = len(df)

    # generate table
    rows = []
    curves_to_plot = ["ExpansionHunter: Q from CIs", "ExpansionHunter: Q from read counts"]
    if not only_expansion_hunter:
        curves_to_plot += ["GangSTR", "HipSTR"]

    for tool in "ExpansionHunter: Q from CIs", "ExpansionHunter: Q from read counts": #, "GangSTR", "HipSTR":
        df_exactly_right = df[df[f"DiffRepeats: Allele: {tool} - Truth"] == 0]

        for min_Q in list(np.arange(0, 0.91, 0.025)) + list(np.arange(0.9, 0.991, 0.01)) + [0.999, 1]:
            total_alleles_passed_Q_threshold = sum(df[f"Q: {tool}"] >= min_Q)
            if total_alleles_passed_Q_threshold == 0:
                continue

            total_alleles_below_Q_threshold = sum(df[f"Q: {tool}"] < min_Q)
            exactly_right_above_Q_threshold = sum(df_exactly_right[f"Q: {tool}"] >= min_Q)
            exactly_right_below_Q_threshold = sum(df_exactly_right[f"Q: {tool}"] < min_Q)

            true_positive_count = exactly_right_above_Q_threshold
            true_negative_count = total_alleles_below_Q_threshold - exactly_right_below_Q_threshold

            rows.append({
                "Q": min_Q,
                "accuracy": (true_positive_count + true_negative_count)/total_allele_count,
                "tool": tool,
                "true positive alleles filtered out": exactly_right_below_Q_threshold/total_allele_count,
                "alleles filtered out": total_alleles_below_Q_threshold/total_allele_count,
            })

    df_plot = pd.DataFrame(rows)
    fig, axes = plt.subplots(2, 1, figsize=(width, height), sharex=True, gridspec_kw={'height_ratios': [1.2, 1.7]})
    fig.suptitle(figure_title, fontsize=TITLE_FONT_SIZE, y=0.97)

    ax0, ax1 = axes
    sns.lineplot(
        data=df_plot,
        x="Q",
        y="true positive alleles filtered out",
        hue="tool",
        ci=None,
        marker="o",
        legend=True,
        ax=ax0,
    )

    sns.lineplot(
        data=df_plot,
        x="Q",
        y="accuracy",
        hue="tool",
        marker="o",
        legend=False,
        ax=ax1,
    )

    ax0.set_xlim(-0.02, 1.02)
    for ax in ax0, ax1:
        ax.set_xticks(np.arange(0.0, 1.01, 0.1))
        ax.xaxis.grid(True, color="#F0F0F0")
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)

    ymax_ax0 = df_plot[df_plot["Q"] < 0.98]["true positive alleles filtered out"].max()
    if not only_expansion_hunter:
        if ymax_ax0 > 0.5:
            ymax_ax0 = 1
        elif ymax_ax0 > 0.25:
            ymax_ax0 = 0.5
        else:
            ymax_ax0 = 0.25
    ax0.set_ylim(0, ymax_ax0 + (0.005 if only_expansion_hunter else 0.02))
    ax0.yaxis.set_major_formatter(lambda y, pos: f"{int(y*total_allele_count):,d} ({100*y:0.0f}%)")
    ax0.set_ylabel("True Positive Alleles Filtered Out", fontsize=13)
    ax0.get_legend().set_title("")
    ax0.get_legend().set_frame_on(False)

    ax1.set_xlabel("Q threshold", fontsize=13)

    ymin_ax1 = df_plot[df_plot["Q"] < (0.99 if only_expansion_hunter else 0.98)].accuracy.min() - 0.01
    if not only_expansion_hunter:
        if ymin_ax1 > 0.75:
            ymin_ax1 = 0.75
        elif ymin_ax1 > 0.5:
            ymin_ax1 = 0.5
        else:
            ymin_ax1 = 0

    ax1.set_ylim(ymin_ax1, 1)
    ax1.set_ylabel("Accuracy", fontsize=13)


def plot_empty_image(figure_title, message, width, height):
    fig, ax = plt.subplots(1, 1, figsize=(width, height))
    fig.suptitle(figure_title, fontsize=TITLE_FONT_SIZE)
    ax.axis('off')
    text = ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=17)

    plt.gcf().canvas.draw()

    bbox = text.get_window_extent()
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    rect = patches.Rectangle(
        bbox.min, bbox.width, bbox.height, linewidth=1, edgecolor='black', facecolor='none', ls='dashed')
    ax.add_patch(rect)


def generate_all_plots(df, args):
    arg_string = " ".join(sys.argv)
    plot_counter = 0

    start_with_plot_i = args.start_with_plot_i
    if start_with_plot_i is None or start_with_plot_i < 0:
        start_with_plot_i = 0

    max_plots = args.n
    if max_plots is None:
        max_plots = 10**9

    output_dir = args.output_dir

    for coverage in ["40x", "30x", "20x", "10x", "exome"] if not args.coverage else [args.coverage]:
        df2 = df[df["Coverage"] == coverage]
        if coverage == "exome":
            df2 = df2[~df2["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

        for motif_size in ["all_motifs", "2bp", "3-6bp", "7-24bp", "25+bp"] if args.min_motif_size is None and args.max_motif_size is None else [f"{args.min_motif_size}-{args.max_motif_size}bp"]:
            if motif_size == "all_motifs":
                df3 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
            elif motif_size == "2bp":
                df3 = df2[df2["MotifSize"] == 2]
            elif motif_size == "3-6bp":
                df3 = df2[(3 <= df2["MotifSize"]) & (df2["MotifSize"] <= 6)]
            elif motif_size == "7-24bp":
                df3 = df2[(7 <= df2["MotifSize"]) & (df2["MotifSize"] <= 24)]
            elif motif_size == "25+bp":
                df3 = df2[(25 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
            else:
                if args.min_motif_size is None and args.max_motif_size is None:
                    raise ValueError(f"Unexpected motif_size value: {motif_size}")

                df3 = df2
                if args.min_motif_size:
                    df3 = df3[(args.min_motif_size <= df3["MotifSize"])]
                if args.max_motif_size:
                    df3 = df3[(df3["MotifSize"] <= args.max_motif_size)]

            for genotype_subset in ["all", "HET", "HOM", "MULTI"] if not args.genotype else [args.genotype]:
                if genotype_subset == "all":
                    df4 = df3
                elif genotype_subset == "HET":
                    df4 = df3[df3["SummaryString"].str.contains(":HET")]
                elif genotype_subset == "HOM":
                    df4 = df3[df3["SummaryString"].str.contains(":HOM")]
                elif genotype_subset == "MULTI":
                    df4 = df3[df3["SummaryString"].str.contains(":MULTI")]
                else:
                    raise ValueError(f"Unexpected genotype_subset value: {genotype_subset}")

                for pure_repeats in ["both", True, False] if "--only-pure-repeats" not in arg_string else [True]:
                    if pure_repeats == "both":
                        df5 = df4
                    elif pure_repeats:
                        df5 = df4[df4["IsPureRepeat"]]
                    else:
                        df5 = df4[~df4["IsPureRepeat"]]

                    hide_no_call_loci_option = "--hide-no-call-loci" in arg_string
                    show_no_call_loci_option = "--show-no-call-loci" in arg_string
                    for exclude_no_call_loci in [False, True] if not hide_no_call_loci_option and not show_no_call_loci_option else ([True] if args.hide_no_call_loci else [False]):
                        if exclude_no_call_loci:
                            if args.only_expansion_hunter:
                                df_plot = df5[~df5["Genotype: ExpansionHunter"].isna()]
                            else:
                                df_plot = df5[
                                    ~df5["Genotype: ExpansionHunter"].isna() &
                                    ~df5["Genotype: GangSTR"].isna() &
                                    ~df5["Genotype: HipSTR"].isna()
                                ]
                        else:
                            df_plot = df5

                        if plot_counter < start_with_plot_i:
                            plot_counter += 1
                            continue

                        print("-"*100)

                        figure_title_line = ""

                        coverage_label = f"exome" if coverage == "exome" else f"{coverage} coverage"

                        filter_description = [f"{coverage_label}"]
                        output_image_filename = "tool_accuracy_vs_Q"
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

                        if exclude_no_call_loci:
                            filter_description.append(f"only loci called by all tools")
                            output_image_filename += f".exclude_no_call_loci"

                        n_locus_ids = len(set(df_plot.LocusId))
                        if len(df_plot) < 50:
                            figure_title_line += f"{n_locus_ids:,d} loci\n\n"+", ".join(filter_description)
                            message = "Not enough alleles to create plot"

                            print(f"Skipping plot {plot_counter}:  {message}")
                            plot_empty_image(figure_title_line, message, args.width, args.height)

                            output_path = os.path.join(output_dir, f"{output_image_filename}.{args.image_type}")
                            plt.savefig(output_path, dpi=300)
                            print(f"Saved {output_path}")
                            plt.close()

                            continue

                        print(f"Generating plot #{plot_counter}")
                        print("Keeping", len(df_plot), "out of ", len(df), "rows")
                        figure_title_line += f"Showing data for {n_locus_ids:,d} loci and {len(df_plot):,d} alleles\n\n" + ", ".join(filter_description) + f""
                        print(figure_title_line)

                        try:
                            plot(df_plot, args.width, args.height,
                                 figure_title=figure_title_line if args.show_title else None,
                                 only_expansion_hunter=args.only_expansion_hunter)

                            output_path = os.path.join(output_dir, f"{output_image_filename}.{args.image_type}")
                            plt.savefig(output_path, bbox_inches='tight', dpi=300)
                            print(f"Saved plot #{plot_counter}: {output_path}")
                            plt.close()

                        except Exception as e:
                            print(f"ERROR: {e}")

                        plot_counter += 1
                        if plot_counter >= start_with_plot_i + max_plots:
                            print(f"Exiting after generating {plot_counter-start_with_plot_i} plot(s)")
                            return


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--start-with-plot-i", type=int, help="If specified, start with this plot number")
    p.add_argument("-n", type=int, help="If specified, only generate this many plots")
    p.add_argument("--width", default=10, type=float)
    p.add_argument("--height", default=10, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"], help="Image type to generate")
    p.add_argument("--show-title", action="store_true", help="Show title in plot")

    g = p.add_argument_group("Filters")
    g.add_argument("--coverage", help="Plot only this coverage (example: \"30x\" or \"exome\")")
    g.add_argument("--min-motif-size", type=int, help="Min motif size")
    g.add_argument("--max-motif-size", type=int, help="Max motif size")
    g.add_argument("--genotype", choices=["all", "HET", "HOM", "MULTI"], help="Plot only this genotype")
    g.add_argument("--only-pure-repeats", action="store_true", help="Plot only loci with pure repeats")
    g2 = g.add_mutually_exclusive_group()
    g2.add_argument("--show-no-call-loci", action="store_true", help="Show loci with no call")
    g2.add_argument("--hide-no-call-loci", action="store_true", help="Hide loci with no call")

    p.add_argument("--only-expansion-hunter", action="store_true", help="Plot only the ExpansionHunter curves")

    p.add_argument("--output-dir", default=".")
    p.add_argument("--verbose", action="store_true", help="Print additional info")
    p.add_argument("combined_tool_results_tsv", nargs="?", default="../tool_comparison/combined.results.alleles.tsv.gz")
    args = p.parse_args()

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)

    if "IsFoundInReference" not in df.columns:
        print("WARNING: IsFoundInReference column not found in input file. Assuming all loci are found in reference...")
        df["IsFoundInReference"] = True

    if "PositiveOrNegative" not in df.columns:
        print("WARNING: PositiveOrNegative column not found in input file. Assuming all loci are positive...")
        df["PositiveOrNegative"] = "positive"

    df = df[df["IsFoundInReference"] & (df["PositiveOrNegative"] == "positive")]

    if args.verbose:
        print("Num loci:")
        print(df.groupby(["PositiveOrNegative", "Coverage", "IsPureRepeat"]).count().LocusId/2)

    print("Computing additional columns...")

    df.loc[:, "Q: ExpansionHunter: Q from CIs"] = 1/np.exp(
        4 * df["CI size: Allele: ExpansionHunter"]/df["NumRepeats: Allele: ExpansionHunter"]
    )

    df.loc[:, "Q: ExpansionHunter: Q from read counts"] = np.where(
        df["FractionOfReadsThatSupportsGenotype: Allele: ExpansionHunter"] < 0.15,
        df["FractionOfReadsThatSupportsGenotype: Allele: ExpansionHunter"] / 0.15,
        1)

    df.loc[:, "DiffRepeats: Allele: ExpansionHunter: Q from CIs - Truth"] = df[
        "DiffRepeats: Allele: ExpansionHunter - Truth"]
    df.loc[:, "DiffRepeats: Allele: ExpansionHunter: Q from read counts - Truth"] = df[
        "DiffRepeats: Allele: ExpansionHunter - Truth"]

    print(f"Generating {str(args.n) + ' ' if args.n else ''}plots",
          f"starting with plot #{args.start_with_plot_i}" if args.start_with_plot_i else "")
    generate_all_plots(df, args)

    print(f"Done")


if __name__ == "__main__":
    main()