import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})


def plot_allele_size_distribution(df, output_dir, plot_type=1, color_by=None, hue_order=None, ax=None, save_image=False):
    # Generate Figure 1: distribution of allele sizes

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 7), dpi=80)

    xlimit = 15 if plot_type == 1 else 39
    ax.xaxis.labelpad = ax.yaxis.labelpad = 15
    ax.set_xlabel(
        "# of Repeats Relative To hg38" if plot_type == 1 else "Allele Size (bp) Relative to hg38",
        fontsize=14,
    )
    ax.set_ylabel("Fraction of Alleles", fontsize=14)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    xticks = range(-xlimit, xlimit + 1, 2) if plot_type == 1 else range(-xlimit, xlimit + 1, 6)
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{x}" if x < 0 else f"+{x}" for x in xticks], fontsize=13)
    ax.set_xlim(-xlimit - 0.52, xlimit + 0.52)

    p = sns.histplot(
        df,
        x="NumRepeatsAwayFromReference" if plot_type == 1 else "NumBasePairsAwayFromReference",
        hue=color_by,
        hue_order=hue_order,
        bins=[b-1.5 for b in range(-xlimit, xlimit + 3, 3)] if plot_type == 2 else None,
        binwidth=None if plot_type == 2 else 1,
        multiple="stack" if not color_by else "fill",
        stat="proportion",
        discrete=True if plot_type == 1 else None,
        ax=ax)

    if plot_type == 1:
        p.set_title("STR Allele Size Distribution", fontsize=14)
        output_image_name = "allele_size_distribution_by_number_of_repeats"
    else:
        p.set_title("STR Allele Size Distribution", fontsize=14)
        output_image_name = "allele_size_distribution_in_base_pairs"

    if color_by:
        output_image_name += f".color_by_{color_by.lower()}"
        if color_by == "Multiallelic":
            ax.get_legend().set_title("Multiallelic Loci")

    plt.yticks(fontsize=13)

    if save_image:
        output_image_path = os.path.join(output_dir, output_image_name)
        ext = ".svg"
        plt.savefig(f"{output_image_path}{ext}", bbox_inches="tight")
        plt.close()
        print(f"Saved {output_image_path}{ext}")

    print(f"Plotted {len(df):,d} allele records")


def plot_motif_distribution(df, output_dir, axes=None, save_image=False):
    print("Plotting allele distribution by motif size")

    hue_limit = 5
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReference"]
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = np.where(
        df["NumRepeatsAwayFromReferenceTruncated"] > hue_limit, hue_limit, df["NumRepeatsAwayFromReferenceTruncated"])
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = np.where(
        df["NumRepeatsAwayFromReferenceTruncated"] < -hue_limit, -hue_limit, df["NumRepeatsAwayFromReferenceTruncated"])
    df = df.sort_values("NumRepeatsAwayFromReferenceTruncated")

    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].astype("int").apply(
        lambda i: f"{i}" if i < 0 else f"+{i}")
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].replace(
        f"-{hue_limit}", f"-{hue_limit} or less")
    df.loc[:, "NumRepeatsAwayFromReferenceTruncated"] = df["NumRepeatsAwayFromReferenceTruncated"].replace(
        f"+{hue_limit}", f"+{hue_limit} or more")

    if axes is None:
        _, axes = plt.subplots(2, 1, figsize=(8, 12), dpi=80)

    ax = axes[0]
    p = sns.histplot(
        df,
        x="NumRepeatsInReference",
        binwidth=1,
        discrete=True,
        ax=ax)

    p.set_title("# of Repeats In hg38 at Truth Set STR Loci", fontsize=14)
    ax.set_xlabel("", fontsize=14)
    ax.set_yscale("log")
    ax.set_ylabel("# of Loci", fontsize=14)

    ax = axes[1]
    sns.histplot(
        df,
        x="NumRepeatsInReference",
        hue="NumRepeatsAwayFromReferenceTruncated",
        palette=sns.color_palette("bwr_r", 10)[::-1],
        binwidth=1,
        discrete=True,
        stat="proportion",
        multiple="fill",
        ax=ax)

    xlimit = 51
    for ax in axes:
        ax.xaxis.labelpad = ax.yaxis.labelpad = 15
        ax.spines.right.set_visible(False)
        ax.spines.top.set_visible(False)
        ax.set_xticks(range(0, xlimit + 1, 5))
        ax.set_xlim(-1, xlimit)

    ax.set_xlabel("# of Repeats in hg38", fontsize=14)
    ax.set_ylabel("Fraction of Alleles", fontsize=14)
    plt.rcParams.update({"text.usetex": True})
    ax.get_legend().set_title("# of Repeats in\nTruth Set Allele\nRelative to hg38\n$_{contractions \quad are < 0}$\n$_{expansions \quad are > 0}   $", prop={'size': 12})
    ax.get_legend()._legend_box.align = "left"
    ax.get_legend().set_bbox_to_anchor((1.35, 0.25))
    ax.get_legend().get_frame().set_color("white")

    if save_image:
        output_image_name = "distribution_by_motif_size"
        output_image_path = os.path.join(output_dir, output_image_name)
        ext = ".svg"
        plt.savefig(f"{output_image_path}{ext}", bbox_inches="tight")
        plt.close()
        print(f"Saved {output_image_path}{ext}")

    print(f"Plotted {len(df):,d} allele records")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--skip-plot1", action="store_true")
    p.add_argument("--skip-plot2", action="store_true")
    p.add_argument("--output-dir", default=".")
    args = p.parse_args()

    input_table_path = "../STR_truthset.v1.alleles.tsv.gz"

    print(f"Reading {input_table_path}")
    df = pd.read_table(input_table_path)
    df = df[df.IsPureRepeat == "Yes"]
    print(f"Parsed {len(df):,d} rows from {input_table_path}")
    df = df.rename(columns={"IsMultiallelic": "Multiallelic"})
    df.loc[:, "NumRepeatsAwayFromReference"] = df.NumRepeats - df.NumRepeatsInReference
    df.loc[:, "NumBasePairsAwayFromReference"] = df.NumRepeatsAwayFromReference * df.MotifSize

    if not args.skip_plot1:
        plot_allele_size_distribution(df, args.output_dir, plot_type=1, save_image=True)
        plot_allele_size_distribution(df, args.output_dir, plot_type=2, save_image=True)
        plot_allele_size_distribution(
            df, args.output_dir, plot_type=1, color_by="Multiallelic", hue_order=["No", "Yes"], save_image=True)

    if not args.skip_plot2:
        plot_motif_distribution(df, args.output_dir, save_image=True)


if __name__ == "__main__":
    main()