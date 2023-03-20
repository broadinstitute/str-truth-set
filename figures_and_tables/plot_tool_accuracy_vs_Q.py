import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import sys

from matplotlib import patches

sns.set_context(font_scale=1.1, rc={
    "font.family": "sans-serif",
})
sns.set_palette(["purple", "blue", "red", "orange"])


FIGURE_SIZE = (10, 10)
TITLE_FONT_SIZE = 14


def plot(df, figure_title=None):
    total_allele_count = len(df)

    # generate table
    rows = []
    for tool in "ExpansionHunter: Q from CIs", "ExpansionHunter: Q from read counts", "GangSTR", "HipSTR":
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
                "alleles filtered out": total_alleles_below_Q_threshold/total_allele_count,
            })

    df_plot = pd.DataFrame(rows)

    fig, axes = plt.subplots(2, 1, figsize=(10, 10), sharex=True, gridspec_kw={'height_ratios': [1.2, 1.7]})
    fig.suptitle(figure_title, fontsize=TITLE_FONT_SIZE, y=0.95)

    ax0, ax1 = axes
    sns.lineplot(
        data=df_plot,
        x="Q",
        y="alleles filtered out",
        hue="tool",
        #style="tool",
        ci=None,
        markers=True,
        legend=True,
        ax=ax0,
    )

    sns.lineplot(
        data=df_plot,
        x="Q",
        y="accuracy",
        hue="tool",
        #style="tool",
        markers=True,
        legend=False,
        ax=ax1,
    )

    ax0.set_xlim(-0.02, 1.02)
    ax0.set_xticks(np.arange(0.0, 1.01, 0.1))
    ax1.set_xticks(np.arange(0.0, 1.01, 0.1))

    ax0.set_ylim(0, df_plot["alleles filtered out"].max())

    ax0.xaxis.grid(True, color="#F0F0F0")
    ax1.xaxis.grid(True, color="#F0F0F0")

    ax0.yaxis.set_major_formatter(lambda y, pos: f"{int(y*total_allele_count):,d} ({100*y:0.0f}%)")

    ax0.set_ylabel(ax0.get_ylabel(), fontsize=13)
    ax0.get_legend().set_title("")
    ax0.get_legend().set_frame_on(False)

    #ax1.set_ylim(df_plot.accuracy.min(), df_plot.accuracy.max())
    ax1.set_ylim(0.5 if df_plot.accuracy.min() > 0.5 else 0, 1)
    ax1.set_xlabel("Q threshold", fontsize=13)
    ax1.set_ylabel(ax1.get_ylabel(), fontsize=13)

    ax0.xaxis.labelpad = ax0.yaxis.labelpad = 15
    ax1.xaxis.labelpad = ax1.yaxis.labelpad = 15
    ax0.spines.right.set_visible(False)
    ax0.spines.top.set_visible(False)
    ax1.spines.right.set_visible(False)
    ax1.spines.top.set_visible(False)


def plot_empty_image(figure_title, message):
    fig, ax = plt.subplots(1, 1, figsize=FIGURE_SIZE, dpi=80)
    fig.suptitle(figure_title, fontsize=TITLE_FONT_SIZE)
    ax.axis('off')
    text = ax.text(0.5, 0.5, message, ha='center', va='center', fontsize=24)

    plt.gcf().canvas.draw()

    bbox = text.get_window_extent()
    bbox = bbox.transformed(fig.dpi_scale_trans.inverted())

    rect = patches.Rectangle(
        bbox.min, bbox.width, bbox.height, linewidth=1, edgecolor='black', facecolor='none', ls='dashed')
    ax.add_patch(rect)


def generate_all_plots(df, output_image_dir, max_plots=None):
    plot_counter = 0
    for motif_size in "all_motifs", "2-6bp", "7-24bp", "25+bp":
        for pure_repeats in "both", True, False:
            for exclude_hipstr_no_call_loci in False, True:
                for coverage in "40x", "30x", "20x", "10x", "exome":
                    for genotype_subset in "all_genotypes", "HET", "HOM", "MULTI":
                        if max_plots is not None and plot_counter >= max_plots:
                            print(f"Exiting after generating {plot_counter} plot(s)")
                            sys.exit(0)

                        print("-"*100)

                        figure_title_line = ""

                        df2 = df.copy()
                        df2 = df2[(df2["Coverage"] == coverage)]

                        if coverage == "exome":
                            df2 = df2[~df2["Genotype: GangSTR"].isna() | ~df2["Genotype: ExpansionHunter"].isna()]
                            df2 = df2[~df2["GeneRegionFromGencode_V42"].isin({"intergenic", "intron", "promoter"})]

                        coverage_label = f"exome" if coverage == "exome" else f"{coverage} coverage"

                        filter_description = [f"{coverage_label}"]
                        output_image_filename = "tool_accuracy_vs_Q"
                        if motif_size == "all_motifs":
                            filter_description.append("all motif sizes")
                            output_image_filename += ".all_motifs"
                            df2 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 50)]
                        elif motif_size == "2-6bp":
                            filter_description.append("2bp to 6bp motifs")
                            output_image_filename += ".2to6bp_motifs"
                            df2 = df2[(2 <= df2["MotifSize"]) & (df2["MotifSize"] <= 6)]
                        elif motif_size == "7-24bp":
                            filter_description.append(f"7bp to 24bp motifs")
                            output_image_filename += ".7to24bp_motifs"
                            df2 = df2[(7 <= df2["MotifSize"]) & (df2["MotifSize"] <= 24)]
                        elif motif_size == "25+bp":
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
                            raise ValueError(f"Unexpected genotype_subset value: {genotype_subset}")

                        if pure_repeats == "both":
                            pass
                        elif pure_repeats:
                            filter_description.append("pure repeats")
                            output_image_filename += ".pure_repeats"
                            df2 = df2[df2["IsPureRepeat"]]
                        else:
                            filter_description.append("repeats with interruptions")
                            output_image_filename += ".with_interruptions"
                            df2 = df2[~df2["IsPureRepeat"]]

                        output_image_filename += f".{coverage}"

                        if exclude_hipstr_no_call_loci:
                            filter_description.append(f"only loci with HipSTR call")
                            output_image_filename += f".exclude_HipSTR_no_call_loci"
                            df2 = df2[~df2["Genotype: HipSTR"].isna()]


                        if len(df2) < 10:
                            message = "Not enough alleles to create plot"

                            figure_title_line += f"{len(set(df2.LocusId)):,d} loci ("+", ".join(filter_description)+")"

                            print(f"Skipping..  {message}")
                            plot_empty_image(figure_title_line, message)
                            plt.savefig(f"{output_image_dir}/{output_image_filename}.png")
                            print(f"Saved {output_image_dir}/{output_image_filename}.png")
                            plt.close()
                            continue

                        print("Keeping", len(df2), "out of ", len(df), "rows")
                        figure_title_line += f"Showing data for {len(set(df2.LocusId)):,d} loci (" + ", ".join(filter_description) + f")"
                        print(figure_title_line)

                        try:
                            plot(df2, figure_title=figure_title_line)

                            plt.savefig(f"{output_image_dir}/{output_image_filename}.png", bbox_inches='tight')
                            print(f"Saved plot #{plot_counter+1}: {output_image_dir}/{output_image_filename}.png")
                            plt.close()

                        except Exception as e:
                            print(f"ERROR: {e}")

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

    df.loc[:, "Q: ExpansionHunter: Q from CIs"] = 1/np.exp(
        4*df["CI size: Allele: ExpansionHunter"]/df["NumRepeats: Allele: ExpansionHunter"]
    )

    df.loc[:, "Q: ExpansionHunter: Q from read counts"] = np.where(
        df["FractionOfReadsThatSupportsGenotype: Allele: ExpansionHunter"] < 0.15,
        df["FractionOfReadsThatSupportsGenotype: Allele: ExpansionHunter"] / 0.15,
        1)

    df.loc[:, "DiffRepeats: Allele: ExpansionHunter: Q from CIs - Truth"] = df[
        "DiffRepeats: Allele: ExpansionHunter - Truth"]
    df.loc[:, "DiffRepeats: Allele: ExpansionHunter: Q from read counts - Truth"] = df[
        "DiffRepeats: Allele: ExpansionHunter - Truth"]

    print("Generating plots...")
    generate_all_plots(df, args.output_dir, max_plots=args.n)

    print(f"Done")


if __name__ == "__main__":
    main()