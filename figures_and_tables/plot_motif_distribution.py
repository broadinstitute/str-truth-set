import argparse
import gzip
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import tqdm

from str_analysis.utils.canonical_repeat_unit import compute_canonical_motif


sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
    "legend.fontsize": 15,
})


def compute_motif_labels(row, existing_motif_labels=(), min_fraction_of_motifs_with_same_length=0.1):
    label = f"{len(row.CanonicalMotif)}bp:{row.CanonicalMotif}"
    if label in existing_motif_labels or (not existing_motif_labels and
        (len(row.CanonicalMotif) <= 2 or
        (row.MotifSize <= 6 and
         row.NumSTRsWithSameCanonicalMotif/row.NumSTRsWithSameMotifSize >= min_fraction_of_motifs_with_same_length))):
        return label
    elif row.MotifSize <= 6:
        return f"{row.MotifSize}bp:other"
    else:
        return f"7+bp"


def plot_hist(df, output_image_filename, title=None, y_label=None, show_title=True, width=10, height=8):
    include_homopolymers = sum(df.MotifSize == 1) > 0

    fig, ax = plt.subplots(figsize=(width, height))
    ax.set_xlim(0, 15.5)
    ax.set_xlabel("Motif Size (bp)", labelpad=15, fontsize=17)
    if y_label:
        ax.set_ylabel(y_label, labelpad=15, fontsize=17)
    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)

    if show_title:
        ax.set_title(title, pad=15, fontsize=17)

    palette = sns.color_palette("hsv", n_colors=len(set(df["MotifLabels"])), desat=0.8)

    sns.histplot(
        df,
        x="MotifSize",
        hue="MotifLabels",
        palette=palette,
        binwidth=1,
        multiple="stack",
        stat="proportion",
        discrete=True,
        legend=True,
        ax=ax)

    # set legend to have 2 columns
    legend = ax.get_legend()
    labels = [text.get_text() for text in legend.get_texts()]
    legend.remove()
    ax.legend(legend.legend_handles, labels, title=None, prop={'size': 16}, frameon=False, ncol=2)

    x_min = 1 if include_homopolymers else 2
    xticks = range(x_min, 16, 1)
    yticks = np.arange(0, 0.75, 0.1)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels([f"{x}" for x in xticks], fontsize=16)
    ax.set_yticklabels([f"{y:0.1f}" for y in yticks], fontsize=16)

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

    p.add_argument("--include-homopolymers", action="store_true", help="Include homopolymers")

    g = p.add_mutually_exclusive_group()
    g.add_argument("--only-pure-repeats", action="store_true", help="Only plot pure repeats")
    g.add_argument("--only-interrupted-repeats", action="store_true", help="Only plot interrupted repeats")

    p.add_argument("--show-title", action="store_true", help="Show title in plot")
    p.add_argument("--width", default=10, type=float)
    p.add_argument("--height", default=8, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"])

    p.add_argument("--truth-set-alleles-table", default="../STR_truth_set.v1.alleles.tsv.gz")
    args = p.parse_args()

    min_ref_bp = "12bp"

    input_bed_path = os.path.join(
        args.ref_dir,
        f"other/repeat_specs_GRCh38_without_mismatches.sorted.trimmed.at_least_{min_ref_bp}.bed.gz",
    )
    rows = []
    print(f"Parsing {input_bed_path}")
    with gzip.open(input_bed_path, "rt") as f:
        for line in tqdm.tqdm(f, unit=" lines", unit_scale=True):
            fields = line.strip().split("\t")
            motif = fields[3]
            if len(motif) == 1 and not args.include_homopolymers:
                continue

            rows.append({
                "chrom": fields[0],
                "start_0based": int(fields[1]),
                "end": int(fields[2]),
                "Motif": motif,
                "num_repeats": float(fields[4]),
                "MotifSize": len(fields[3]),
                "CanonicalMotif": compute_canonical_motif(fields[3], include_reverse_complement=True),
            })

    print(f"Parsed {len(rows):,d} rows")
    df_ref = pd.DataFrame(rows)
    print("Adding columns...")
    df_ref = add_columns(df_ref)

    print(f"Plotting motifs from {len(df_ref):,d} reference loci...")
    plot_hist(df_ref, os.path.join(args.output_dir, f"motif_distribution_in_hg38_with_atleast_{min_ref_bp}.{args.image_type}"),
              y_label="Fraction of All Pure TR Loci in hg38",
              title=f"Motif Distribution in hg38" if args.show_title else None,
              width=args.width,
              height=args.height,
    )

    image_name = "motif_distribution"
    df = pd.read_table(args.truth_set_alleles_table)
    if not args.include_homopolymers:
        df = df[df["MotifSize"] > 1]
    else:
        image_name += ".including_homopolymers"

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
              y_label="Fraction of Truth Set Alleles",
              title="Motif Distribution in Truth Set" if args.show_title else None,
              width=args.width,
              height=args.height,
    )


if __name__ == "__main__":
    main()