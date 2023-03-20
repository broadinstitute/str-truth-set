import collections
import pandas as pd
from figures_and_tables.numbers_utils import format_np

df_variants = pd.read_table("pure_STR_truth_set.v1.variants.tsv.gz")
df_variants = df_variants[df_variants.IsPureRepeat]

#print("\n".join(df_variants.columns))
columns = ["OverlapsIlluminaSTRCatalog: Locus",
           "OverlapsGangSTRCatalog13: Locus", "OverlapsGangSTRCatalog17: Locus",
           "OverlapsHipSTRCatalog: Locus"] + [f"OverlapsTRFPureRepeats{i}bp: Locus" for i in (6, 9, 12, 15)]


print(format_np(len(df_variants[df_variants.MotifSize <= 6]), len(df_variants)), " of variants have 2-6bp motifs")
for label, df in [
    ("all", df_variants), ("2-6bp", df_variants[df_variants.MotifSize <= 6])
]:
    print("="*100)

    print("    -----")
    print(f"% STR truth set loci ({label} motifs) missed by different catalogs: ")
    for column in columns:
        count_values = collections.Counter(df[column])
        print(format_np(count_values["no overlap"], len(df)), f"of STR variants ({label} motifs) don't overlap", column.replace("Overlaps", "").split(":")[0])

    

    print("    -----")
    print(f"% STR truth set loci ({label} motifs) missed by different catalogs for loci with a reference locus >= 24bp: ")

    df_filtered = df[df.NumRepeatsInReference * df.MotifSize >= 24]
    print(format_np(len(df_filtered), len(df)), "have a reference locus size of at least 24bp")
    for column in columns:
        count_values = collections.Counter(df_filtered[column])
        print(format_np(count_values["no overlap"], len(df_filtered)), f"of STR variants ({label} motifs) don't overlap", column.replace("Overlaps", "").split(":")[0])


    

    print("    -----")
    print(f"% STR truth set loci ({label} motifs) missed by different catalogs for loci with a reference locus >= 5 repeats: ")

    df_filtered = df[(df.NumRepeatsLongAllele - df.NumRepeatsInReference) >= 5]
    print(format_np(len(df_filtered), len(df)), "have a reference locus >= 5 repeats")
    for column in columns:
        count_values = collections.Counter(df_filtered[column])
        print(format_np(count_values["no overlap"], len(df_filtered)), f"of STR variants ({label} motifs) don't overlap", column.replace("Overlaps", "").split(":")[0])

    

    print("    -----")
    print(f"% STR truth set loci ({label} motifs) missed by different catalogs for loci with a reference locus >= 10 repeats: ")

    df_filtered = df[(df.NumRepeatsLongAllele - df.NumRepeatsInReference) >= 10]
    print(format_np(len(df_filtered), len(df)), "have a reference locus >= 10 repeats")
    for column in columns:
        count_values = collections.Counter(df_filtered[column])
        print(format_np(count_values["no overlap"], len(df_filtered)), f"of STR variants ({label} motifs) don't overlap", column.replace("Overlaps", "").split(":")[0])
