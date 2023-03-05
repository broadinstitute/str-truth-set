import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import numpy as np
from matplotlib.ticker import MultipleLocator

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})


def generate_plot(with_coverage=False):

    if with_coverage:
        input_path = "../tool_comparison/STR_tool_timing.with_coverage.tsv"
    else:
        input_path = "../tool_comparison/STR_tool_timing.tsv"

    df = pd.read_table(input_path)

    df = df[df["positive_or_negative_loci"] == "positive_loci"]

    if with_coverage:
        df = df[~df["tool"].isin({"IlluminaExpansionHunter",})]
        df = df[~df["coverage"].isin({"exome", "5x",})]

    df.loc[:, "tool"] = df.tool.replace({
        "ExpansionHunter":         "(Optimized) ExpansionHunter v5",
        "IlluminaExpansionHunter": "(Original) ExpansionHunter v5",
        "GangSTR":                 "GangSTR v2.5",
        "HipSTR":                  "HipSTR v0.6.2",
    })

    if with_coverage:
        df.loc[:, "coverage"] = df.coverage.replace({"40x genome": "40x"})
        df.loc[:, "tool_and_coverage"] = df["tool"] + " (" + df["coverage"] + ")"

        sort_by = "tool_and_coverage"
        order_dict = {
            f"{tool_label} ({coverage})": (i, j) for i, tool_label in enumerate(
                ["(Original) ExpansionHunter v5", "(Optimized) ExpansionHunter v5", "GangSTR v2.5", "HipSTR v0.6.2"]
            ) for j, coverage in enumerate([
                "40x", "30x", "20x", "10x", "5x", "exome",
            ])
        }
    else:
        sort_by = "tool"
        order_dict = {
            v: i for i, v in enumerate(
                ["(Original) ExpansionHunter v5", "(Optimized) ExpansionHunter v5", "GangSTR v2.5", "HipSTR v0.6.2"]
            )
        }

    df = df.sort_values(by=sort_by, key=lambda tool_column: [order_dict[t] for t in tool_column])

    df.loc[:, "max_resident_set_kbytes"] = df["max_resident_set_kbytes"] / 10**6

    if with_coverage:
        figsize = (22, 10)
        x_column = "tool_and_coverage"
        title_suffix = "\nby genome coverage"
        y_limit1 = 100
        y_ticks1 = np.arange(0, 101, 10)
        output_image_name_suffix = ".by_coverage"
    else:
        figsize = (12, 7)
        x_column = "tool"
        title_suffix = ""
        y_limit1 = 300
        y_ticks1 = np.arange(0, 301, 25)
        output_image_name_suffix = ""

    # generate plot
    plt.rc('font', family='Verdana')
    fig, axes = plt.subplots(figsize=figsize, ncols=2, tight_layout=True)

    for i, ax in enumerate(axes):
        if i == 0:
            title = f"Run Time{title_suffix}"
            y_column = "minutes_per_10k_loci"
            y_label = "Minutes per 10,000 loci"
            y_limit = y_limit1
            y_ticks = y_ticks1
            y_minor_ticks = 25
            y_major_ticks = 50
        elif i == 1:
            title = f"Memory Usage{title_suffix}"
            y_column = "max_resident_set_kbytes"
            y_label = "Max memory used (Gb)"
            y_limit = 2
            y_ticks = np.arange(0, 2, 0.25)
            y_minor_ticks = 0.125
            y_major_ticks = 0.25
        else:
            raise ValueError(f"i: {i}")

        sns.stripplot(
            x=x_column,
            y=y_column,
            data=df,
            jitter=True,
            color="#555555",
            ax=ax)

        # plot the mean line
        boxplot = sns.boxplot(
            x=x_column,
            y=y_column,
            data=df,
            showmeans=False,
            meanline=False,
            medianprops={'visible': True, 'color': 'red', 'ls': '-', 'lw': 2},
            whiskerprops={'visible': False},
            zorder=10,
            showfliers=False,
            showbox=False,
            showcaps=False,
            width=0.4,
            ax=ax)

        ax.spines.top.set_visible(False)
        ax.spines.right.set_visible(False)

        ax.set_title(title, fontsize=16, pad=20, fontweight=500)
        ax.set_xlabel(None)
        ax.set_ylabel(y_label, labelpad=15, fontsize=14)

        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=14)

        ax.set_ylim(0, y_limit)
        ax.yaxis.set_ticks(y_ticks)
        ax.yaxis.set_major_locator(MultipleLocator(y_major_ticks))
        ax.yaxis.set_minor_locator(MultipleLocator(y_minor_ticks))

        ax.yaxis.grid(True)
        ax.yaxis.grid(which='minor', color='lightgray', linewidth=0.5, linestyle=':', alpha=0.9)
        ax.yaxis.grid(which='major', color='lightgray', linewidth=0.5, linestyle='-', alpha=0.9)

        if with_coverage:
            ax.vlines([3.5, 7.5], ymin=0, ymax=ax.get_ylim()[-1])

        # label medians
        medians = df.groupby([x_column])[y_column].median()

        lines = boxplot.axes.get_lines()
        for xtick in boxplot.axes.get_xticks():
            y = lines[2+xtick*3].get_ydata()[0]
            print(xtick, medians[xtick])
            boxplot.axes.text(
                xtick + (0.5 if i == 0 else 0.6), y, int(y) if i == 0 else round(y, 2),
                horizontalalignment='center',
                verticalalignment='center',
                color='black',
                fontsize=12,
                )
        print("---")

    print(f"Plotted {len(df):,d} records")

    output_image_name = f"tool_runtime_and_memory{output_image_name_suffix}.svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    print(f"Saved {output_image_name}")
    plt.close()


def main():
    generate_plot(with_coverage=False)
    generate_plot(with_coverage=True)


if __name__ == "__main__":
    main()


