import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

sns.set_context(font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "legend.title_fontsize": 14,
    "axes.labelsize": 13,
    "axes.titlesize": 15,
})


def plot_mutation_rate_by_allele_size(df, args, show_homopolymers=False, show_dinucleotides=False):
    """This function plots the fraction of alleles that are multi-allelic (a proxy for mutation rates) with allele size
    on the x-axis, and color based on the motif size.

    df (pd.DataFrame):
    """
    if show_homopolymers:
        df = df[(df["MotifSize"] >= 1) & (df["MotifSize"] <= 3) & (df["CanonicalMotif"] != "CG")]
    elif show_dinucleotides:
        df = df[(df["MotifSize"] >= 2) & (df["MotifSize"] <= 3) & (df["CanonicalMotif"] != "CG")]
    else:
        df = df[(df["MotifSize"] >= 3) & (df["MotifSize"] <= 6)]

    x_column = "RepeatSize (bp)"
    x_column_binned = "Total Allele Size (bp)"
    bin_size = 10

    df = df[df.FractionPureRepeats >= 0.75]

    if show_homopolymers:
        df["hue"] = df["MotifSize"].astype(str) + "bp" + np.where(
            df["MotifSize"] <= 2, " (" + df["CanonicalMotif"] + "), ", ", ") + np.where(
            df["IsPureRepeat"], "pure", "interrupted")
        df = df[df["IsPureRepeat"]]
    elif show_dinucleotides:
        df["hue"] = df["MotifSize"].astype(str) + "bp" + np.where(
            df["MotifSize"] == 2, " (" + df["CanonicalMotif"] + "), ", ", ") + np.where(
            df["IsPureRepeat"], "pure", "interrupted")
    else:
        df["hue"] = df["MotifSize"].astype(str) + "bp, " + np.where(df["IsPureRepeat"], "pure", "interrupted")

    df[x_column_binned] = pd.cut(df[x_column], list(range(0, 100, bin_size))).apply(lambda x: x.mid)
    df["DataPoints"] = 1

    df2 = df.groupby([x_column_binned, "hue", "MotifSize", "IsPureRepeat"]).agg(
        {"IsMultiallelic": "mean", "DataPoints": "sum"}
    ).reset_index()

    df2 = df2[(df2["DataPoints"] >= 20) & (df2[x_column_binned].astype(int) < 65)]
    df2.sort_values(by=["IsPureRepeat", "MotifSize"], ascending=[False, True], inplace=True)

    if show_homopolymers:
        hue_order = [
            "1bp (C), pure",
            "1bp (A), pure",
            "2bp (AT), pure",
            "2bp (AC), pure",
            "2bp (AG), pure",
            "3bp, pure",
            #"2bp (AT), interrupted",
            #"2bp (AC), interrupted",
            #"2bp (AG), interrupted",
            #"3bp, interrupted",
        ]

        palette = list(sns.color_palette("Purples_d", 2))[::-1]
        palette += list(sns.color_palette("Oranges", 3))[::-1] + ["#52EC68"]
        palette += list(sns.color_palette("Blues", 3))[::-1] + ["#c8f9cf"]
    elif show_dinucleotides:
        hue_order = [
            "2bp (AT), pure",
            "2bp (AC), pure",
            "2bp (AG), pure",
            "3bp, pure",
            "2bp (AT), interrupted",
            "2bp (AC), interrupted",
            "2bp (AG), interrupted",
            "3bp, interrupted",
        ]
        "52EC68"  # 3bp
        "53BDEC"  # 4bp
        "9D53EC"  # 5bp
        "EC53A8"  # 6bp


        palette = list(sns.color_palette("Oranges", 3))[::-1] + ["#52EC68"]
        palette += list(sns.color_palette("Blues", 3))[::-1] + ["#c8f9cf"]
    else:
        hue_order = None
        #palette = list(sns.color_palette("Oranges", len(set(df2[df2["IsPureRepeat"]]["MotifSize"]))))
        #palette += list(sns.color_palette("Blues", len(set(df2[~df2["IsPureRepeat"]]["MotifSize"]))))
        palette = [
            "#52EC68",  # 3bp
            "#53BDEC",  # 4bp
            "#9D53EC",  # 5bp
            "#EC53A8",  # 6bp
            "#c8f9cf",  # lightened by 68%
            "#c8eaf9",
            "#e0c8f9",
            "#f9c8e3",
        ]

    fig, ax = plt.subplots(figsize=(args.width, args.height))
    sns.lineplot(data=df2, x=x_column_binned, y="IsMultiallelic", hue="hue", hue_order=hue_order,
                 palette=palette, marker="o", legend="full", ax=ax)
    ax.set_xlabel(x_column_binned, labelpad=15)
    ax.set_ylabel("Fraction Multi-allelic Variants", labelpad=15)
    ax.set_ylim(-0.05, 0.82)

    ax.get_legend().set_title(None)
    ax.get_legend().set_frame_on(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    # set padding between axis label and tick labels
    ax.tick_params(axis="both", which="major", pad=10)

    # save figure
    filename = "mutation_rates_by_allele_size"
    if show_homopolymers:
        filename += ".1bp_and_2bp_motifs"
    elif show_dinucleotides:
        filename += ".2bp_motifs"
    else:
        filename += ".3-6bp_motifs"
    filename += f".{args.image_type}"

    output_path = os.path.join(args.output_dir, filename)
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    print(f"Saved {output_path}")


def plot_mutation_rate_by_fraction_interrupted_repeats(df, args):
    #df = df[(df["MotifSize"] >= 2) & (df["MotifSize"] <= 6)]
    df = df[(df["MotifSize"] >= 3) & (df["MotifSize"] <= 6)]

    x_column = "FractionInterruptedRepeats"
    x_column_binned = "Fraction Interrupted Repeats"
    bin_size = 0.1

    df[x_column_binned] = pd.cut(df[x_column], list(np.arange(0, 1, bin_size))).apply(lambda x: x.mid)
    #df["hue"] = np.where(df["MotifSize"] == 2, "2bp", "3-6bp")

    df["DataPoints"] = 1
    #df2 = df.groupby([x_column_binned]).agg(
    df2 = df.groupby([x_column_binned, "MotifSize"]).agg(
        {"IsMultiallelic": "mean", "DataPoints": "sum"}
    ).reset_index()
    df2 = df2[(df2["DataPoints"] >= 20)]

    fig, ax = plt.subplots(figsize=(args.width/2 + 2, args.height))
    sns.lineplot(data=df2, x=x_column_binned, y="IsMultiallelic", marker="o", ax=ax)  # hue="hue", legend="full",

    #ax.set_xlim(0, 0.75)
    ax.set_xlim(0, 0.5)
    ax.set_xlabel("Fraction Interrupted Repeats", labelpad=15)
    ax.set_ylabel("Fraction Multi-allelic Variants", labelpad=15)
    #ax.get_legend().set_title("Motif Size")
    #ax.get_legend().set_frame_on(False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # save figure
    output_path = os.path.join(args.output_dir, f"mutation_rates_by_fraction_interrupted_repeats.{args.image_type}")
    fig.savefig(output_path, bbox_inches="tight", dpi=300)
    print(f"Saved {output_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default=".")
    p.add_argument("--image-type", default="svg", choices=["svg", "png"], help="Image type to save")
    p.add_argument("--width", default=12, type=float, help="Width of image")
    p.add_argument("--height", default=8, type=float, help="Height of image")
    p.add_argument("--truth-set-alleles-table", default="../STR_truth_set.v1.alleles.tsv.gz")
    args = p.parse_args()

    df = pd.read_table(args.truth_set_alleles_table)
    df["NumInterruptedRepeats"] = df.NumRepeats - df.NumPureRepeats
    df["FractionInterruptedRepeats"] = df.NumInterruptedRepeats/df.NumRepeats

    print(f"Plotting mutation rates based on {len(df):,d} truth set alleles")
    plot_mutation_rate_by_allele_size(df, args, show_dinucleotides=False)
    plot_mutation_rate_by_allele_size(df, args, show_dinucleotides=True)
    plot_mutation_rate_by_allele_size(df, args, show_homopolymers=True)

    plot_mutation_rate_by_fraction_interrupted_repeats(df, args)

    mu_pure = df[df.IsPureRepeat].IsMultiallelic.mean()
    print(f"Mutation rate for pure repeats: {mu_pure:.3f}")
    mu_interrupted = df[~df.IsPureRepeat].IsMultiallelic.mean()
    print(f"Mutation rate for interrupted repeats: {mu_interrupted:.3f}")
    print(f"Ratio of mutation rates: {mu_pure/mu_interrupted:.3f}")


if __name__ == "__main__":
    main()
