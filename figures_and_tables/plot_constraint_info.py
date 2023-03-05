import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})


def plot_gene_contraint(df):
    fig, axes = plt.subplots(1, 3, sharex=False, sharey=True, figsize=(20, 6), dpi=80)

    for column, ax in zip([
        'pLI',
        "O/E LoF upperbound (LOEUF)",
        "O/E missense upperbound",
        #'pRecessive',
    ], axes):

        sns.violinplot(
            y="Source",
            x=column,
            data=df,
            inner=None,
            cut=0,
            scale="width",
            ax=ax)

        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        for violin in ax.collections:
            bbox = violin.get_paths()[0].get_extents()
            x0, y0, width, height = bbox.bounds
            violin.set_clip_path(plt.Rectangle((x0, y0), width, height / 2, transform=ax.transData))

        old_len_collections = len(ax.collections)
        stripplot_data = df[~df.Source.isin(["All Genes"])]
        stripplot_data = stripplot_data.sort_values('InheritanceMode & AgeOfOnset')
        sns.stripplot(
            x=column,
            y="Source",
            data=stripplot_data,
            hue='InheritanceMode & AgeOfOnset',
            ax=ax,
            palette=["gray", "red", "orange", "blue", "blue"],
        )

        for dots in ax.collections[old_len_collections:]:
            dots.set_offsets(dots.get_offsets() + np.array([0, 0.2]))
        if column == "pLI":
            i = 0
            for _, row in stripplot_data[~stripplot_data.gene_name.isna()].iterrows():
                text_x = row.pLI
                text_y = 1.05
                if row.gene_name == "ATXN3":
                    text_x += 0.02
                    text_y += 0.15
                elif row.gene_name == "PABPN1":
                    text_x -= 0.08
                    text_y += 0.35
                elif row.gene_name == "SOX3":
                    text_x += 0.02
                    text_y += 0.08
                elif row.gene_name == "HOXD13":
                    text_x -= 0.02
                    text_y += 0.4
                else:
                    continue
                ax.text(x=text_x, y=text_y, s=row.gene_name)

        sns.boxplot(
            y="Source", x=column, data=df,
            saturation=1, showfliers=False, width=0.1,
            medianprops={'zorder': 5, 'linewidth': 3},
            whiskerprops={'zorder': 5, 'linewidth': 3},
            boxprops={'zorder': 3, 'linewidth': 3, 'facecolor': 'white'},
            ax=ax)

        if column == "pLI":
            ax.get_legend().set_title(f"")
            sns.move_legend(ax, loc=(-1, 0))
        else:
            ax.get_legend().remove()

        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel("")

    output_image_name = "gene_constraint.svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    print(f"Saved {output_image_name}")
    plt.close()


def main():
    constraint_df = pd.read_table("../ref/other/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz", compression="gzip")
    constraint_df = constraint_df[["gene_id", "pLI", "oe_lof_upper", "oe_mis_upper", "oe_syn_upper", "pRec", "mu_syn"]]
    constraint_df = constraint_df.set_index("gene_id")

    known_pathogenic_df = pd.read_table("../ref/other/private/known_pathogenic_STRs_with_Utah_columns.tsv")
    known_pathogenic_df = known_pathogenic_df[~known_pathogenic_df["GeneId"].isna()]
    known_pathogenic_df = known_pathogenic_df.set_index("GeneId")
    known_pathogenic_df = known_pathogenic_df[known_pathogenic_df["Gene Region"].str.startswith("coding:")]
    known_pathogenic_df = known_pathogenic_df[known_pathogenic_df["RU length"] <= 6]

    truth_set_df = pd.read_table("../STR_truth_set.v1.variants.tsv.gz")
    truth_set_df = truth_set_df[truth_set_df["IsPureRepeat"] == "Yes"]
    truth_set_df = truth_set_df.set_index("GeneIdFromMane_V1")
    truth_set_df = truth_set_df[truth_set_df.GeneRegionFromMane_V1 == "CDS"]
    print(f"Loaded {len(truth_set_df)} truth set loci that MANE v1 CDS")

    df1 = constraint_df[constraint_df.index.isin(known_pathogenic_df.index)].reset_index()

    df1 = df1.set_index("gene_id").join(known_pathogenic_df[
                                            ["InheritanceMode", 'age_onset_min: Utah', 'age_onset_max: Utah', "Gene"]
                                        ], how="left")
    df1 = df1.reset_index().rename(columns={
        "index": "gene_id",
        "Gene": "gene_name",
        "age_onset_min: Utah": "age_onset_min",
        "age_onset_max: Utah": "age_onset_max",
    })

    def combine_inhertiance_and_age_of_onset(row):
        if row.InheritanceMode == "AD" or  row.InheritanceMode == "XD":
            return row.InheritanceMode + (", late onset" if (row["age_onset_max"] + row["age_onset_min"])/2 > 25 else "")
        else:
            return row.InheritanceMode

    df1.loc[:, "InheritanceMode & AgeOfOnset"] = df1.apply(combine_inhertiance_and_age_of_onset, axis=1)

    df2 = constraint_df[constraint_df.index.isin(truth_set_df.index)].reset_index()
    print(len(df2), "truth set STRs in MANE v1 CDSes have constraint info")
    df3 = constraint_df

    df1.loc[:, "Source"] = "Genes Containing\nDisease Associated STR Loci"
    df2.loc[:, "Source"] = "Genes Containing\nTruth Set STR Variants"
    df3.loc[:, "Source"] = "All Genes"
    df2.loc[:, "InheritanceMode"] = df2.loc[:, "InheritanceMode & AgeOfOnset"] = ""

    df_final = pd.concat([df1, df2, df3])

    df_final = df_final.sort_values("Source", ascending=False)
    df_final = df_final.rename(columns={
        "pRec": "pRecessive",
        "oe_lof_upper": "O/E LoF upperbound (LOEUF)",
        "oe_mis_upper": "O/E missense upperbound",
    })

    plot_gene_contraint(df_final)


if __name__ == "__main__":
    main()