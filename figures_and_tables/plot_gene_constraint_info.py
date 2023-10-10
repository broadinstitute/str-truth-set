import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
from numbers_utils import format_np

np.random.seed(10)

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "axes.labelsize": 15,
})


DISEASE_ASSOCIATED_STRs_LABEL = "Genes Containing\nDisease Associated TR Loci"
TRUTH_SET_TRs_LABEL = "Genes Containing\nTruth Set TR Alleles"
ALL_GENES_LABEL = "All Genes"

Y_AXIS_ORDER = {
    DISEASE_ASSOCIATED_STRs_LABEL: 3,
    TRUTH_SET_TRs_LABEL: 2,
    ALL_GENES_LABEL: 1,
}

def plot_gene_constraint(df, args):
    fig, axes = plt.subplots(1, 3, sharex=False, sharey=True, figsize=(20, 6))
    if args.show_title:
        suptitle_artist = fig.suptitle("Gene Constraint Metrics and TR Alleles", fontsize=16, y=0.97)
        extra_artists = [suptitle_artist]
    else:
        extra_artists = []

    df.index = range(len(df))
    y_axis_order = [key for key in sorted(Y_AXIS_ORDER, key=Y_AXIS_ORDER.get)]
    for column_i, (x_column, ax) in enumerate(zip([
        'pLI',
        "O/E LoF upperbound (LOEUF)",
        "O/E missense upperbound",
        #'pRecessive',
    ], axes)):
        ax.xaxis.labelpad = 15
        sns.violinplot(
            y="Source",
            x=x_column,
            data=df,
            inner=None,
            cut=0,
            # select colors
            palette=["#3A9239", "#3274A1", "#E1812B"],
            scale="width",
            order=y_axis_order,
            ax=ax)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        for violin in ax.collections:
            bbox = violin.get_paths()[0].get_extents()
            x0, y0, width, height = bbox.bounds
            violin.set_clip_path(plt.Rectangle((x0, y0 - 0.1), width, height / 2 + 0.1, transform=ax.transData))

        for strip_plot_category in TRUTH_SET_TRs_LABEL, DISEASE_ASSOCIATED_STRs_LABEL:
            old_len_collections = len(ax.collections)
            stripplot_data = df[df.Source.isin([strip_plot_category])]
            stripplot_data = stripplot_data.sort_values('InheritanceMode')
            sns.stripplot(
                x=x_column,
                y="Source",
                data=stripplot_data,
                hue='InheritanceMode',
                order=y_axis_order,
                ax=ax,
                palette=["gray"] if strip_plot_category == TRUTH_SET_TRs_LABEL else ["red", "blue", "orange"],
                alpha=0.5,
            )

            for dots in ax.collections[old_len_collections:]:
                dots.set_offsets(dots.get_offsets() + np.array([0, 0.2]))

            if strip_plot_category == DISEASE_ASSOCIATED_STRs_LABEL:
                for _, row in stripplot_data[~stripplot_data.gene_name.isna()].iterrows():
                    if row.gene_name not in {"RUNX2", "HOXD13"}:
                        continue

                    text_x = row[x_column] + 0.02
                    text_y = 2.4
                    ax.text(x=text_x, y=text_y, s=row.gene_name)

        sns.boxplot(
            y="Source", x=x_column, data=df,
            saturation=1, showfliers=False, width=0.1,
            medianprops={'zorder': 5, 'linewidth': 3},
            whiskerprops={'zorder': 5, 'linewidth': 3},
            boxprops={'zorder': 3, 'linewidth': 3, 'facecolor': 'white'},
            ax=ax)

        if column_i == 1:
            ax.get_legend().set_title(f"")
            sns.move_legend(ax, loc=(0.55, 0.16))
        else:
            ax.get_legend().remove()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel("")

    output_path = os.path.join(args.output_dir, f"gene_constraint_metrics_and_STRs.{args.image_type}")
    plt.savefig(output_path, bbox_extra_artists=extra_artists, bbox_inches="tight", dpi=300)

    print(f"Saved {output_path}")
    plt.close()


def plot_str_variant_size_vs_gene_constraint(df, args):
    fig, axes = plt.subplots(1, 3, sharex=False, sharey=True, figsize=(20, 6))
    if args.show_title:
        suptitle_artist = fig.suptitle("Gene Constraint Metrics by TR Allele Size", fontsize=16, y=0.97)
        extra_artists = [suptitle_artist]
    else:
        extra_artists = []

    for column, ax in zip([
        'pLI',
        "O/E LoF upperbound (LOEUF)",
        "O/E missense upperbound",
        #'pRecessive',
    ], axes):
        sns.scatterplot(data=df, x=column, y="TR Allele Size (bp)", ax=ax)
        xlim = list(ax.get_xlim())
        xlim[0] = -0.05
        ax.set_xlim(xlim)
        ax.set_ylim([-100, 100])
        ax.grid(which='major', linewidth=0.5)

    output_path = os.path.join(args.output_dir, f"gene_constraint_metrics_vs_STR_allele_size.{args.image_type}")
    plt.savefig(output_path, bbox_extra_artists=extra_artists, bbox_inches="tight", dpi=300)

    print(f"Saved {output_path}")
    plt.close()



def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default=".")
    p.add_argument("--show-title", action="store_true")
    p.add_argument("--image-type", default="svg", choices=["svg", "png"])

    p.add_argument("--constraint-table", default="../ref/other/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
    p.add_argument("--disease-associated-loci-table", default="../ref/other/known_disease_associated_STR_loci.tsv")
    p.add_argument("--truth-set-alleles-table", default="../STR_truth_set.v1.alleles.tsv.gz")
    args = p.parse_args()

    constraint_df = pd.read_table(args.constraint_table, compression="gzip")
    constraint_df = constraint_df[["gene_id", "pLI", "oe_lof_upper", "oe_mis_upper", "oe_syn_upper", "pRec", "mu_syn"]]
    constraint_df = constraint_df.rename(columns={
        "pRec": "pRecessive",
        "oe_lof_upper": "O/E LoF upperbound (LOEUF)",
        "oe_mis_upper": "O/E missense upperbound",
    })
    constraint_df = constraint_df.set_index("gene_id")
    print(f"Loaded {len(constraint_df):,d} rows from the gnomAD gene constraint table")

    known_pathogenic_df = pd.read_table(args.disease_associated_loci_table)
    known_pathogenic_df = known_pathogenic_df[~known_pathogenic_df["GeneId"].isna()]
    known_pathogenic_df = known_pathogenic_df[known_pathogenic_df["Gene Region"].str.startswith("coding:")]
    known_pathogenic_df = known_pathogenic_df[known_pathogenic_df["Motif Length"] <= 6]
    known_pathogenic_df = known_pathogenic_df.set_index("GeneId")
    known_pathogenic_df.loc[:, "InheritanceMode"] = known_pathogenic_df["InheritanceMode"].replace({
        "AD": "Autosomal Dominant",
        "AR": "Autosomal Recessive",
        "XR": "X-linked Recessive",
    })
    print(f"Loaded {len(known_pathogenic_df):,d} coding known disease-associated STR loci")

    truth_set_df = pd.read_table(args.truth_set_alleles_table)
    truth_set_df = truth_set_df[truth_set_df["IsPureRepeat"]]
    truth_set_df = truth_set_df.set_index("GeneIdFromMane_V1")
    truth_set_df = truth_set_df[truth_set_df.GeneRegionFromMane_V1 == "CDS"]
    print(f"Loaded {len(truth_set_df):,d} truth set loci that overlap MANE v1 coding regions")


    print("----")
    df1 = constraint_df[constraint_df.index.isin(known_pathogenic_df.index)]
    print(f"{len(df1):,d} genes with constraint info contain coding known disease-associated STR loci")

    df1 = df1.join(known_pathogenic_df[["InheritanceMode", "Gene"]], how="left")

    df1 = df1.reset_index().rename(columns={
        "index": "gene_id",
        "Gene": "gene_name",
    })

    df2 = constraint_df[constraint_df.index.isin(truth_set_df.index)].reset_index()
    print(f"{len(df2):,d} genes with constraint info contain coding truth set STR alleles")

    df3 = constraint_df.reset_index()
    print(f"{len(df3):,d} genes with constraint info total")

    df1.loc[:, "Source"] = DISEASE_ASSOCIATED_STRs_LABEL
    df2.loc[:, "Source"] = TRUTH_SET_TRs_LABEL
    df3.loc[:, "Source"] = ALL_GENES_LABEL

    df2.loc[:, "InheritanceMode"] = ""

    df_final = pd.concat([df1, df2, df3])
    df_final.loc[:, "SourceOrder"] = df_final["Source"].map(Y_AXIS_ORDER)
    df_final = df_final.sort_values("SourceOrder", ascending=True)

    print("----")
    print("pLI constraint score max:", df_final.pLI.max())
    print("LOEUF constraint score max:", df_final["O/E LoF upperbound (LOEUF)"].max())
    print("missense constraint score max:", df_final["O/E missense upperbound"].max())
    print("----")

    print(f"{format_np(sum(df3.pLI > 0.9), len(df3))} all genes have pLI > 0.9")
    print(f"{format_np(sum(df2.pLI > 0.9), len(df2))} ({sum(df2.pLI > 0.9)/len(df2)/(sum(df3.pLI > 0.9)/len(df3)):0.1f} above background) truth set variant genes have pLI > 0.9")
    print(f"{format_np(sum(df1.pLI > 0.9), len(df1))} ({sum(df1.pLI > 0.9)/len(df1)/(sum(df3.pLI > 0.9)/len(df3)):0.1f} above background) known disease-associated genes have pLI > 0.9")

    print(f"{format_np(sum(df3['O/E LoF upperbound (LOEUF)'] < 0.35), len(df3))} all genes have LEOUF < 0.35")
    print(f"{format_np(sum(df2['O/E LoF upperbound (LOEUF)'] < 0.35), len(df2))} ({sum(df2['O/E LoF upperbound (LOEUF)'] < 0.35)/len(df2)/(sum(df3['O/E LoF upperbound (LOEUF)'] < 0.35)/len(df3)):0.1f} above background) truth set variant genes have LEOUF < 0.35")
    print(f"{format_np(sum(df1['O/E LoF upperbound (LOEUF)'] < 0.35), len(df1))} ({sum(df1['O/E LoF upperbound (LOEUF)'] < 0.35)/len(df1)/(sum(df3['O/E LoF upperbound (LOEUF)'] < 0.35)/len(df3)):0.1f} above background) known disease-associated genes have LEOUF < 0.35")

    plot_gene_constraint(df_final, args)

    truth_set_df = truth_set_df[truth_set_df.GeneRegionFromMane_V1 == "CDS"]
    truth_set_df.loc[:, 'NumRepeatsDiffFromReference'] = truth_set_df.NumRepeats - truth_set_df.NumRepeatsInReference
    truth_set_df.loc[:, 'TR Allele Size (bp)'] = truth_set_df.NumRepeatsDiffFromReference * truth_set_df.MotifSize
    truthset_with_constraint_df = truth_set_df.join(constraint_df, how="inner")

    truthset_with_constraint_df = truthset_with_constraint_df.reset_index()
    truthset_with_constraint_df.index = range(len(truthset_with_constraint_df))
    plot_str_variant_size_vs_gene_constraint(truthset_with_constraint_df, args)

    print("---")
    constraint_df = pd.read_table(args.constraint_table, compression="gzip")
    print(f"Genes with pLI > 0.9 are {constraint_df[constraint_df.pLI > 0.9].cds_length.median()/constraint_df[constraint_df.pLI < 0.9].cds_length.median():0.1f}x longer than genes with pLI < 0.9 (based on median)")
    print(f"Genes with pLI > 0.9 are {constraint_df[constraint_df.pLI > 0.9].cds_length.mean()/constraint_df[constraint_df.pLI < 0.9].cds_length.mean():0.1f}x longer than genes with pLI < 0.9 (based on mean)")

if __name__ == "__main__":
    main()