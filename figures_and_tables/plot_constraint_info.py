import gzip
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from tqdm import tqdm

sns.set_context("paper", font_scale=1.1, rc={
    "font.family": "sans-serif",
})


def plot_distribution(df):

    fig, ax = plt.subplots(figsize=(20, 3.5), dpi=80, tight_layout=True)

    print(f"Plotted {len(df):,d} allele records")

    output_image_name = "syndip_indel_size_distribution.svg"
    plt.savefig(f"{output_image_name}", bbox_inches="tight")
    print(f"Saved {output_image_name}")
    plt.close()


def main():
    plot_distribution(df)


if __name__ == "__main__":
    main()