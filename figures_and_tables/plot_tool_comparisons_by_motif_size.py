import argparse
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})

GREEN_COLOR = "#50AA44"


def plot_distribution_by_motif_size(df, figure_title, output_image_path, show_title=True):
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
            ax.set_xlim(left=min(df_current["MotifSize"])-0.5, right=max(df_current["MotifSize"])+0.5)

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
                x="MotifSize",
                hue=hue_column,
                hue_order=["ExactlyTheSame", "OverlappingCIs", "Discordant"],
                binwidth=1,
                palette=[GREEN_COLOR, "orange", "#FF665533"], # "#5588FF33"
                multiple="stack" if row_i == 0 else "fill",
                stat="proportion",
                discrete=True,
                legend=(column_j == 0) and (row_i > 0),
                ax=axes[column_j])

            if show_title:
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
    if show_title:
        suptitle_artist = fig.suptitle(figure_title, fontsize=20, y=1.02)
        extra_artists = [suptitle_artist]
    else:
        extra_artists = []

    plt.savefig(f"{output_image_path}", bbox_extra_artists=extra_artists, bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_path}")


def generate_all_distribution_by_motif_size_plots(df, output_image_dir, max_plots=None):
    plot_counter = 0
    for coverage in "40x", "30x", "20x", "10x", "05x", "exome":
        df_current = df[df["coverage"] == coverage]

        coverage_label = f"Exome Data" if coverage == "exome" else f"{coverage} Coverage WGS Data"

        output_image_filename = "tool_accuracy_by_motif_size"
        output_image_filename += ".pure_repeats"
        output_image_filename += f".{coverage}"

        output_image_path = f"{output_image_dir}/{output_image_filename}.svg"
        df_current = df_current.sort_values("RepeatSize (bp): Allele: Truth")
        plot_distribution_by_motif_size(df_current, f"Accuracy by Motif Size ({coverage_label})", output_image_path)

        plot_counter += 1
        if max_plots is not None and plot_counter >= max_plots:
            print(f"Exiting after generating {plot_counter} plot(s)")
            return


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", type=int, help="If specified, only generate this many plots. Useful for testing")
    p.add_argument("--output-dir", default="../tool_comparison/figures")
    p.add_argument("combined_tool_results_tsv", nargs="?", default="../tool_comparison/combined.results.alleles.tsv.gz")
    args = p.parse_args()

    print(f"Loading {args.combined_tool_results_tsv}")
    df = pd.read_table(args.combined_tool_results_tsv)
    df = df[df["IsFoundInReference"] &
            (df["PositiveOrNegative"] == "positive") &
            (df["MotifSize"] >= 2) & (df["MotifSize"] <= 24)]

    print("Generating plots...")
    generate_all_distribution_by_motif_size_plots(df, args.output_dir, max_plots=args.n)

    print(f"Done")


if __name__ == "__main__":
    main()