import argparse
import os
import pandas as pd


CONCORDANCE_PAIRS = [
    ("ExpansionHunter", "Truth"),
    ("GangSTR", "Truth"),
    ("ExpansionHunter", "GangSTR"),
]

warning_counter = 0

def compute_concordance_label_func_wrapper(column_suffix1="ExpansionHunter", column_suffix2="GangSTR"):
    """Creates function for adding concordance columns"""

    def compute_concordance_between_calls(row):
        global warning_counter

        label1 = column_suffix1
        label2 = column_suffix2
        concordance_labels = []

        label_map = {
            0: "Discordant",
            1: "OverlappingCIs",
            2: "ExactlyTheSame",
        }

        for allele_number in 1, 2:
            if pd.isna(row[f"NumRepeats: Allele {allele_number}: {label1}"]) and pd.isna(row[f"NumRepeats: Allele {allele_number}: {label2}"]):
                #raise ValueError("Both tools have NumRepeats == NaN")
                warning_counter += 1
                if not row[f"IsHomRef: {label1}"] and not row[f"IsHomRef: {label2}"]:
                    print(f"WARNING {warning_counter}: Both {label1} and {label2} have NumRepeats: Allele {allele_number}: {label1} and Allele {allele_number}: {label2} == NaN")
                continue

            elif pd.isna(row[f"NumRepeats: Allele {allele_number}: {label1}"]):
                if row[f"IsHomRef: {label2}"]:
                    concordance_labels.append(2)
                else:
                    concordance_labels.append(0)
            elif pd.isna(row[f"NumRepeats: Allele {allele_number}: {label2}"]):
                if row[f"IsHomRef: {label1}"]:
                    concordance_labels.append(2)
                else:
                    concordance_labels.append(0)
            elif row[f"NumRepeats: Allele {allele_number}: {label1}"] == row[f"NumRepeats: Allele {allele_number}: {label2}"]:
                concordance_labels.append(2)
            elif row[f"CI end: Allele {allele_number}: {label1}"] >= row[f"CI start: Allele {allele_number}: {label2}"] and \
                row[f"CI start: Allele {allele_number}: {label1}"] <= row[f"CI end: Allele {allele_number}: {label2}"]:
                # |--i1---|
                #     |----i2---| i2_end
                concordance_labels.append(1)
            else:
                concordance_labels.append(0)

        if len(concordance_labels) == 0:
            return label_map[2], label_map[2], label_map[2]

        concordance_label = min(concordance_labels)

        # return short allele concordance, long allele concordance, and combined concordance
        return label_map[concordance_labels[0]], label_map[concordance_labels[-1]], label_map[concordance_label]

    return compute_concordance_between_calls


# function for adding concordance column
def compute_max_distance_func_wrapper(column_suffix1="ExpansionHunter", column_suffix2="GangSTR"):

    def compute_max_distance_between_calls(row):
        diffs = []
        for allele_number in 1, 2:
            if pd.isna(row[f"NumRepeats: Allele {allele_number}: {column_suffix1}"]):
                #print(f"NumRepeats: Allele {allele_number}: {column_suffix1} is NA")
                continue
            if pd.isna(row[f"NumRepeats: Allele {allele_number}: {column_suffix2}"]):
                #print(f"NumRepeats: Allele {allele_number}: {column_suffix2} is NA")
                continue
            num_repeats_tool1_minus_tool2 = int(row[f"NumRepeats: Allele {allele_number}: {column_suffix1}"]) - int(row[f"NumRepeats: Allele {allele_number}: {column_suffix2}"])
            diffs.append(num_repeats_tool1_minus_tool2)

        if not diffs:
            return None

        if abs(diffs[0]) > abs(diffs[1]):
            return diffs[0]
        else:
            return diffs[1]

    return compute_max_distance_between_calls


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--debug", action="store_true", help="Whether to print additional info about input and output columns.")
    p.add_argument("--output-tsv", help="Output path of combined tsv file")
    p.add_argument("combined_tsv", help="Path of the combined tsv containing GangSTR, EH and other results.")
    
    args = p.parse_args()

    if not os.path.isfile(args.combined_tsv):
        p.error(f"{args.combined_tsv} not found")

    return args


def main():
    args = parse_args()

    df = pd.read_table(args.combined_tsv)

    print("Adding concordance columns...")
    for label1, label2 in CONCORDANCE_PAIRS:
        concordance_column_name = f"Concordance: {label1} vs {label2}"
    
        df[[f"Allele 1: {concordance_column_name}", f"Allele 2: {concordance_column_name}", concordance_column_name]] = \
            df.apply(compute_concordance_label_func_wrapper(label1, label2), axis=1, result_type="expand")
    
        print(f"Computed {concordance_column_name} column. {sum(df[concordance_column_name].isna())} of {len(df)} values are NA")

    print("Adding max diff columns...")
    for label1, label2 in CONCORDANCE_PAIRS:
        max_diff_column_name = f"Max Diff: {label1} - {label2}"
        df.loc[:, max_diff_column_name] = df.apply(compute_max_distance_func_wrapper(label1, label2), axis=1)
    
        for allele_number in 1, 2:
            diff_column_name = f"Diff: Allele {allele_number}: {label1} - {label2}"
            #df.loc[:, max_diff_column_name] = int(row[f"NumRepeats: Allele {allele_number}: {column_suffix1}"]) - int(row[f"NumRepeats: Allele {allele_number}: {column_suffix2}"])
    
        print(f"Computed {max_diff_column_name} column. {sum(df[concordance_column_name].isna())} of {len(df)} values are NA")
    
    if "Coverage: Truth" in set(df.columns):
        df = df.drop("Coverage: Truth", axis=1)

    if not args.output_tsv:
        args.output_tsv = args.combined_tsv.replace(".tsv", ".with_concordance.tsv")

    df.to_csv(args.output_tsv, index=False, header=True, sep="\t")

    print("Output columns:")
    for column in sorted(df.columns):
        is_nan_count = sum(pd.isna(df[column]))
        print(f"\t{is_nan_count:7,d} of {len(df):,d} ({100*is_nan_count/len(df):5.1f}%) NaN values in {column}")

    print(f"Wrote {len(df):,d} rows to {args.output_tsv}")


if __name__ == "__main__":
    main()
