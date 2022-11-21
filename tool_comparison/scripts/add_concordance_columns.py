import argparse
import os
import pandas as pd


CONCORDANCE_PAIRS = [
    ("ExpansionHunter", "Truth"),
    ("GangSTR", "Truth"),
    ("ExpansionHunter", "GangSTR"),
]

WARNING_COUNTER = 0

def compute_concordance_label_func_wrapper(column_suffix1="ExpansionHunter", column_suffix2="GangSTR"):
    """Creates function for adding concordance columns"""

    def compute_concordance_between_calls(row):
        global WARNING_COUNTER

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
                WARNING_COUNTER += 1
                if not row[f"IsHomRef: {label1}"] and not row[f"IsHomRef: {label2}"]:
                    print(f"WARNING {WARNING_COUNTER:3d}: Both {label1} and {label2} have NumRepeats: Allele {allele_number}: {label1} and Allele {allele_number}: {label2} == NaN")
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
def compute_distance_func_wrapper(column_suffix1="ExpansionHunter", column_suffix2="GangSTR"):
    def convert_num_repeats_to_string(num_repeats):
        if num_repeats > 0:
            sign = "+"
        elif num_repeats < 0:
            sign = "-"
        else:
            sign = ""

        if abs(num_repeats) < 5:
            return f"{sign}{num_repeats}"
        elif abs(num_repeats) <= 9:
            return f"{sign}5 to {sign}9"
        elif abs(num_repeats) <= 14:
            return f"{sign}10 to {sign}14"
        elif abs(num_repeats) <= 19:
            return f"{sign}15 to {sign}19"
        elif abs(num_repeats) <= 24:
            return f"{sign}20 to {sign}24"
        elif abs(num_repeats) <= 50:
            return f"{sign}25 to {sign}49"
        elif abs(num_repeats) <= 100:
            return f"{sign}50 to {sign}99"
        else:
            return f"{sign}100 or more"

    def convert_base_pairs_to_string(base_pairs):
        if base_pairs > 0:
            sign = "+"
        elif base_pairs < 0:
            sign = "-"
        else:
            sign = ""

        if abs(base_pairs) <= 12:
            return f"{sign}1 to {sign}12bp"
        elif abs(base_pairs) <= 24:
            return f"{sign}13 to {sign}24bp"
        elif abs(base_pairs) <= 48:
            return f"{sign}25 to {sign}48bp"
        elif abs(base_pairs) <= 72:
            return f"{sign}48 to {sign}72bp"
        elif abs(base_pairs) <= 96:
            return f"{sign}72 to {sign}96bp"
        elif abs(base_pairs) <= 144:
            return f"{sign}96 to {sign}144bp"
        elif abs(base_pairs) <= 192:
            return f"{sign}144 to {sign}192bp"
        else:
            return f"{sign}193bp or more"

    def compute_distance_between_calls(row):
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

        diff1 = diffs[0]
        diff2 = diffs[1]

        return (
            [convert_num_repeats_to_string(d) for d in (diff1, diff2)] +
            [convert_base_pairs_to_string(d * row.MotifSize) for d in (diff1, diff2)]
        )

    return compute_distance_between_calls


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
        if "truth" not in label1.lower() and "truth" not in label2.lower():
            continue
        #max_diff_column_name = f"Max Diff: {label1} - {label2}"
        diff_repeats1_column_name = f"DiffRepeats: Allele 1: {label1} - {label2}"
        diff_repeats2_column_name = f"DiffRepeats: Allele 2: {label1} - {label2}"
        diff_size1_column_name = f"DiffSize: Allele 1: {label1} - {label2}"
        diff_size2_column_name = f"DiffSize: Allele 2: {label1} - {label2}"
        df[[diff_repeats1_column_name, diff_repeats2_column_name, diff_size1_column_name, diff_size2_column_name]] = \
            df.apply(compute_distance_func_wrapper(label1, label2), axis=1, result_type="expand")

        for column_name in diff_repeats1_column_name, diff_repeats2_column_name, diff_size1_column_name, diff_size2_column_name:
            print(f"Computed {column_name} column. {sum(df[column_name].isna())} of {len(df)} values are NA")

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
