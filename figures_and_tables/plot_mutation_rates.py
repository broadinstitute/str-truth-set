import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns


def plot(df):
    """This function plots the fraction of alleles that are multi-allelic (a proxy for mutation rates) with allele size
    on the x-axis, and color based on the motif size.

    df (pd.DataFrame):
    """


def add_columns(df):
    df.loc[:, "IsMultiallelic"]
    return df


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default=".")
    p.add_argument("--ref-dir", default="../ref")
    p.add_argument("--truth-set-variants-table", default="../STR_truth_set.v1.variants.tsv.gz")
    args = p.parse_args()

    df = pd.read_table(args.truth_set_variants_table)
    df = add_columns(df)
    print(f"Plotting mutation rates from {len(df):,d} truth set loci")


#if __name__ == "__main__":
#    main()

#%%

df = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df.columns

#%%

df["MotifSize"]  # integer
df["IsMultiallelic"] # boolean
df["IsPureRepeat"] # boolean
df["RepeatSize (bp)"]  # integer
df["NumPureRepeats"]  # integer
df["NumRepeats"]  # integer
#df["RepeatSizeShortAllele (bp)"] # integer
#df["RepeatSizeLongAllele (bp)"]  # integer


# this code uses seaborn to create a line plot with the x-axis being the repeat size, and the y-axis being the fraction of alleles that are multi-allelic
# the color of the line is based on the motif size

#%%

#
#x_column = "RepeatSizeLongAllele (bp)"
x_column, x_column_binned, bin_size = "RepeatSize (bp)", "Allele Size (bp)", 5
#x_column, x_column_binned, bin_size = "NumRepeats", "# of Repeats", 4

df_filtered = df[(df["MotifSize"] >= 3) & (df["MotifSize"] <= 6) & (df[x_column] > 0) & (df[x_column] < 100)]

df_filtered[x_column_binned] = pd.cut(df_filtered[x_column], list(range(0, 100, bin_size))).apply(lambda x: x.mid)
df_filtered["DataPoints"] = 1
df2 = df_filtered.groupby([x_column_binned, "MotifSize", "IsPureRepeat"]).agg({"IsMultiallelic": "mean", "DataPoints": "sum"}).reset_index()
df2["hue"] = df2["MotifSize"].astype(str) + "bp " + np.where(df2["IsPureRepeat"], "pure", "impure")

df2 = df2[df2["DataPoints"] >= 100]
df2.sort_values(by=["IsPureRepeat", "MotifSize"], ascending=[False, True], inplace=True)

# palette with 5 colors that are green (darkest to lightest)

palette = list(sns.color_palette("Oranges", len(set(df2[df2["IsPureRepeat"]]["MotifSize"]))))
# palette with 4 colors that are blue (darkest to lightest)
palette += list(sns.color_palette("Blues", len(set(df2[~df2["IsPureRepeat"]]["MotifSize"]))))

# scatter plot connected with lines
ax = sns.lineplot(data=df2, x=x_column_binned, y="IsMultiallelic", hue="hue", palette=palette, marker="o", legend="full")
ax.set(xlabel=x_column_binned, ylabel="Fraction of alleles that are multi-allelic",
       title="Mutation Rates in the Truth Set")
ax.get_legend().set_title(None)
plt.show()

# NumPureRepeats
df2

#%%
#%%

#TODO use most common motif for interrupted repeats

