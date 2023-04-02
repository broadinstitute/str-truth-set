import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
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


def plot_hist(df, output_image_filename, title=None, y_label=None, show_title=True, width=10, height=8):
    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlim(0, 15.5)
    ax.set_xlabel("Motif size (bp)", labelpad=15, fontsize=15)
    if y_label:
        ax.set_ylabel(y_label, labelpad=15, fontsize=15)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    if show_title:
        ax.set_title(title, pad=15, fontsize=15)

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
    ax.set_xticklabels([f"{x}" for x in xticks], fontsize=14)
    ax.set_yticklabels([f"{y:0.1f}" for y in yticks], fontsize=14)

    plt.savefig(f"{output_image_filename}", bbox_inches="tight", dpi=300)
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
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default=".")
    p.add_argument("--ref-dir", default="../ref")

    g = p.add_mutually_exclusive_group()
    g.add_argument("--only-pure-repeats", action="store_true", help="Only plot pure repeats")
    g.add_argument("--only-interrupted-repeats", action="store_true", help="Only plot interrupted repeats")

    p.add_argument("--show-title", action="store_true", help="Show title in plot")
    p.add_argument("--width", default=10, type=float)
    p.add_argument("--height", default=8, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"])

    p.add_argument("--truth-set-variants-table", default="../STR_truth_set.v1.variants.tsv.gz")
    args = p.parse_args()

    min_ref_bp = "12bp"
    df_ref = pd.read_table(os.path.join(args.ref_dir,
                           f"other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_{min_ref_bp}.bed.gz"),
        names=["chrom", "start_0based", "end", "Motif", "num_repeats"],
        index_col=False)
    df_ref.loc[:, "MotifSize"] = df_ref.Motif.apply(lambda m: len(m))
    df_ref.loc[:, "CanonicalMotif"] = df_ref.Motif.apply(lambda m: compute_canonical_motif(m, include_reverse_complement=True))
    df_ref = add_columns(df_ref)

    print(f"Plotting motifs from {len(df_ref):,d} reference loci")
    plot_hist(df_ref, os.path.join(args.output_dir, f"motif_distribution_in_hg38_with_atleast_{min_ref_bp}.{args.image_type}"),
              y_label="Fraction of TR Loci in hg38",
              title=f"Motif Distribution in hg38" if args.show_title else None,
              width=args.width,
              height=args.height,
    )

    image_name = "motif_distribution"
    df = pd.read_table(args.truth_set_variants_table)
    if args.only_pure_repeats:
        df = df[df["IsPureRepeat"]]
        print(f"Plotting motifs from {len(df):,d} pure truth set loci")
        image_name += ".pure_repeats"
    elif args.only_interrupted_repeats:
        df = df[~df["IsPureRepeat"]]
        print(f"Plotting motifs from {len(df):,d} truth set loci with interruptions")
        image_name += ".interrupted_repeats"

    df = add_columns(df, existing_motif_labels=set(df_ref.MotifLabels))
    print(f"Plotting motifs from {len(df):,d} pure truth set loci")
    plot_hist(df, os.path.join(args.output_dir, f"{image_name}.{args.image_type}"),
              y_label="Fraction of Truth Set Loci",
              title="Motif Distribution in Truth Set" if args.show_title else None,
              width=args.width,
              height=args.height,
    )


if __name__ == "__main__":
    main()