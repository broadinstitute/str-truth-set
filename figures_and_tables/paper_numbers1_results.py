
import collections
import pandas as pd
from figures_and_tables.numbers_utils import format_n, format_np, search

df_variants = pd.read_table("STR_truth_set.v1.variants.tsv.gz")
df_variants_with_interruptions = df_variants[~df_variants.IsPureRepeat]
df_variants = df_variants[df_variants.IsPureRepeat]

df_alleles = pd.read_table("STR_truth_set.v1.alleles.tsv.gz")
df_alleles_with_interruptions = df_alleles[~df_alleles.IsPureRepeat]
df_alleles = df_alleles[df_alleles.IsPureRepeat]
