import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "legend.fontsize": 12,
})


def compute_motif_labels(row, existing_motif_labels=(), min_fraction_of_motifs_with_same_length=0.1):
    label = f"{len(row.CanonicalMotif)}bp:{row.CanonicalMotif}"
    if label in existing_motif_labels or (not existing_motif_labels and
        (len(row.CanonicalMotif) == 2 or
        (row.MotifSize <= 6 and
         row.NumSTRsWithSameCanonicalMotif/row.NumSTRsWithSameMotifSize >= min_fraction_of_motifs_with_same_length))):
        return label
    elif row.MotifSize <= 6:
        return f"{row.MotifSize}bp:other"
    else:
        return f"7+bp"


def plot_hist(df, output_image_filename, title=None, y_label=None, show_title=True):
    fig, ax = plt.subplots(figsize=(8, 9), dpi=80)
    ax.set_xlim(0, 15.5)
    ax.set_xlabel("Motif size (bp)", labelpad=15, fontsize=14)
    if y_label:
        ax.set_ylabel(y_label, labelpad=15, fontsize=14)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    if show_title:
        ax.set_title(title, pad=15, fontsize=14)

    palette = sns.color_palette("hsv", n_colors=len(set(df["MotifLabels"])), desat=0.8)

    df = df.rename(columns={
        "MotifLabels": "Normalized Motifs",
    })

    sns.histplot(
        df,
        x="MotifSize",
        hue="Normalized Motifs",
        palette=palette,
        binwidth=1,
        multiple="stack",
        stat="proportion",
        discrete=True,
        ax=ax)

    ax.get_legend().set_title("Normalized Motif", prop={'size': 13})
    ax.get_legend().set_frame_on(False)

    xticks = range(2, 16, 1)
    yticks = np.arange(0, 0.75, 0.1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{x}" for x in xticks], fontsize=13)
    ax.set_yticklabels([f"{y:0.1f}" for y in yticks], fontsize=13)

    plt.savefig(f"{output_image_filename}", bbox_inches="tight")
    plt.close()
    print(f"Saved {output_image_filename}")


def add_columns(df, existing_motif_labels=()):
    df.loc[:, "NumSTRsWithSameCanonicalMotif"] = df.groupby('CanonicalMotif')['CanonicalMotif'].transform('count')
    df.loc[:, "NumSTRsWithSameMotifSize"] = df.groupby('MotifSize')['MotifSize'].transform('count')
    df.loc[:, "MotifLabels"] = df.apply(lambda r: compute_motif_labels(r, existing_motif_labels=existing_motif_labels, min_fraction_of_motifs_with_same_length=0.07), axis=1)
    df.loc[:, "MoreMotifLabels"] = df.apply(lambda r: compute_motif_labels(r, existing_motif_labels=existing_motif_labels, min_fraction_of_motifs_with_same_length=0.05), axis=1)
    df.sort_values(["MotifSize", "NumSTRsWithSameCanonicalMotif", "Motif"], ascending=[True, False, True], inplace=True)
    return df


def main():
    min_ref_bp = "12bp"
    df_ref = pd.read_table(
        f"../ref/other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_{min_ref_bp}.bed.gz",
        names=["chrom", "start_0based", "end", "Motif", "num_repeats"],
        index_col=False)
    df_ref.loc[:, "MotifSize"] = df_ref.Motif.apply(lambda m: len(m))
    df_ref.loc[:, "CanonicalMotif"] = df_ref.Motif.apply(lambda m: compute_canonical_motif(m, include_reverse_complement=True))
    df_ref = add_columns(df_ref)
    print(f"Plotting motifs from {len(df_ref):,d} reference loci")
    plot_hist(df_ref, f"motif_distribution_in_hg38_with_atleast_{min_ref_bp}.svg",
              y_label="Fraction of Pure STR Loci in hg38",
              title=f"Motifs of Pure STR Loci in hg38")

    df = pd.read_table("../STR_truth_set.v1.variants.tsv.gz")
    df_pure = df[df["IsPureRepeat"] == "Yes"]
    df_pure = add_columns(df_pure, existing_motif_labels=set(df_ref.MotifLabels))
    print(f"Plotting motifs from {len(df_pure):,d} pure truth set loci")
    plot_hist(df_pure, f"motif_distribution.pure_repeats.svg",
              y_label="Fraction of Pure STR Variants in Truth Set",
              title="Motifs of STR Variants In Truth Set")

    df_with_interruptions = df[df["IsPureRepeat"] == "No"]
    df_with_interruptions = add_columns(df_with_interruptions)
    print(f"Plotting motifs from {len(df_with_interruptions):,d} truth set loci with interruptions")
    plot_hist(df_with_interruptions, f"motif_distribution.with_interruptions.svg",
              y_label="Fraction of Interrupted STR Variants in Truth Set",
              title="Motifs of STR Variants In Truth Set\n(only interrupted repeats)")


if __name__ == "__main__":
    main()