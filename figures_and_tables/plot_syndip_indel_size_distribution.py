import argparse
import gzip
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from tqdm import tqdm

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
    "svg.fonttype": "none",  # add text as text rather than curves
})


def plot_distribution(df, args):
    print(f"Rendering indel size distribution with width: {args.width}")
    fig, ax = plt.subplots(figsize=(args.width, args.height), dpi=80, tight_layout=True)

    ax.spines.right.set_visible(False)
    ax.spines.top.set_visible(False)
    ax.set_xlabel("CHM1-CHM13 Deletion & Insertion Size (Kb)", labelpad=15, fontsize=16)
    ax.set_ylabel("Allele Count", labelpad=15, fontsize=16)
    ax.set_yscale('log')
    ax.set_xlim((-20, 20))
    sns.histplot(
        df,
        x="indel_kb",
        binwidth=0.1,
        color="#0033FF",
        ax=ax)

    fig.tight_layout()

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=14)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=15)
    ax.xaxis.set_major_formatter('{x:0.0f} Kb')

    print(f"Plotted {len(df):,d} allele records")

    output_path = os.path.join(args.output_dir, f"syndip_indel_size_distribution.{args.image_type}")
    plt.savefig(output_path, bbox_inches="tight", dpi=300)
    print(f"Saved {output_path}")
    plt.close()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--syndip-vcf", default="../step1.high_confidence_regions.vcf.gz")

    p.add_argument("--width", default=20, type=float)
    p.add_argument("--height", default=3.5, type=float)
    p.add_argument("--image-type", default="svg", choices=["svg", "png"])
    p.add_argument("--output-dir", default=".")
    args = p.parse_args()

    indels = set()
    indel_sizes = []

    with gzip.open(args.syndip_vcf, "rt") as f:
        for line in tqdm(f, unit=" lines", total=4_081_580):
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alts = fields[4].split(",")
            for alt in alts:
                indel_size = len(alt) - len(ref)
                if indel_size == 0:
                    continue
                if (chrom, pos, ref, alt) in indels:
                    continue
                indels.add((chrom, pos, ref, alt))
                indel_sizes.append(indel_size)

    df = pd.DataFrame({"indel_size": indel_sizes})
    df.loc[:, "indel_kb"] = df.indel_size / 10**3
    df.indel_size.mean()

    print("Largest insertion size: ", max(df.indel_kb), "kb")
    print("Largest deletion size: ", min(df.indel_kb), "kb")
    print("Mean event size: ", df.indel_size.mean(), "kb")

    plot_distribution(df, args)


if __name__ == "__main__":
    main()