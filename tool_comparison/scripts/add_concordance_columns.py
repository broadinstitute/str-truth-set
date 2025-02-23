import argparse
import os
import pandas as pd

#CONCORDANCE_PAIRS = [
    #("ExpansionHunter", "Truth"),
    #("GangSTR", "Truth"),
    #("ExpansionHunter", "GangSTR"),
    #("HipSTR", "Truth"),
    #("TRGT", "Truth"),
#]

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
            elif all(c in row for c in (
                f"CI start: Allele {allele_number}: {label1}", f"CI end: Allele {allele_number}: {label1}",
                f"CI start: Allele {allele_number}: {label2}", f"CI end: Allele {allele_number}: {label2}",
            )) and row[f"CI end: Allele {allele_number}: {label1}"] >= row[f"CI start: Allele {allele_number}: {label2}"] and \
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
    """Creates function for computing distance between two genotypes"""

    def compute_distance_between_calls(row):
        diffs = []
        for allele_number in 1, 2:
            if pd.isna(row[f"NumRepeats: Allele {allele_number}: {column_suffix1}"]) or \
               pd.isna(row[f"NumRepeats: Allele {allele_number}: {column_suffix2}"]):
                return None, None, None, None

            num_repeats_tool1_minus_tool2 = int(row[f"NumRepeats: Allele {allele_number}: {column_suffix1}"]) - \
                                            int(row[f"NumRepeats: Allele {allele_number}: {column_suffix2}"])
            diffs.append(num_repeats_tool1_minus_tool2)

        return diffs[0], diffs[1], diffs[0] * row.MotifSize, diffs[1] * row.MotifSize

    return compute_distance_between_calls


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--debug", action="store_true", help="Whether to print additional info about input and output columns.")
    p.add_argument("--output-tsv", help="Output path of combined tsv file")
    p.add_argument("--tool", help="Which tool to compare to the true genotype", required=True,
        choices={"ExpansionHunter", "ExpansionHunter-dev", "GangSTR", "HipSTR", "TRGT", "LongTR", "NewTruthSet"})
    p.add_argument("--compare-to", help="Which tool to compare to the true genotype", default="Truth",
        choices={"Truth", "ExpansionHunter", "ExpansionHunter-dev", "GangSTR", "HipSTR", "LongTR", "TRGT"})
    p.add_argument("combined_tsv", help="Path of the combined tsv containing GangSTR, EH and other results.")
    
    args = p.parse_args()

    if not os.path.isfile(args.combined_tsv):
        p.error(f"{args.combined_tsv} not found")

    return args


def main():
    args = parse_args()

    df = pd.read_table(args.combined_tsv)

    df.rename(columns={
        "LocusSize (bp)": "ReferenceLocusSize (bp)",
    }, inplace=True)

    label1 = args.tool
    label2 = args.compare_to

    print("Adding concordance columns...")
    concordance_column_name = f"Concordance: {label1} vs {label2}"

    df[[f"Allele 1: {concordance_column_name}", f"Allele 2: {concordance_column_name}", f"Variant: {concordance_column_name}"]] = \
        df.apply(compute_concordance_label_func_wrapper(label1, label2), axis=1, result_type="expand")

    concordance_column_name = f"Variant: {concordance_column_name}"
    print(f"Computed {concordance_column_name} column. {sum(df[concordance_column_name].isna())} of {len(df)} values are NA")

    print("Adding max diff columns...")
    #max_diff_column_name = f"Max Diff: {label1} - {label2}"
    diff_repeats1_column_name = f"DiffRepeats: Allele 1: {label1} - {label2}"
    diff_repeats2_column_name = f"DiffRepeats: Allele 2: {label1} - {label2}"
    diff_size1_column_name = f"DiffSize (bp): Allele 1: {label1} - {label2}"
    diff_size2_column_name = f"DiffSize (bp): Allele 2: {label1} - {label2}"
    df[[diff_repeats1_column_name, diff_repeats2_column_name, diff_size1_column_name, diff_size2_column_name]] = \
        df.apply(compute_distance_func_wrapper(label1, label2), axis=1, result_type="expand")

    for column_name in diff_repeats1_column_name, diff_repeats2_column_name, diff_size1_column_name, diff_size2_column_name:
        print(f"Computed {column_name} column. {sum(df[column_name].isna())} of {len(df)} values are NA")

    if not args.output_tsv:
        args.output_tsv = args.combined_tsv.replace(".tsv", ".with_concordance.tsv")

    df = df.astype(str).replace(r"\.0$", "", regex=True)
    df.to_csv(args.output_tsv, index=False, header=True, sep="\t")

    print("Output columns:")
    for column in sorted(df.columns):
        is_nan_count = sum(pd.isna(df[column]))
        print(f"\t{is_nan_count:7,d} of {len(df):,d} ({100*is_nan_count/len(df):5.1f}%) NaN values in {column}")

    print(f"Wrote {len(df):,d} rows to {args.output_tsv}")

    write_alleles_table(df, args.output_tsv.replace(".tsv", ".alleles.tsv"))


def write_alleles_table(df, output_tsv_path):
    output_columns = []
    for key in df.columns:
        if "Allele 1" in key and "Allele 1: " not in key:
            # assuming that all allele-specific columns include substring "Allele n: "
            raise ValueError(f"Unexpected column name: {key}")
        if "Allele 2" in key and "Allele 2: " not in key:
            # assuming that all allele-specific columns include substring "Allele n: "
            raise ValueError(f"Unexpected column name: {key}")

        if "Allele 1: " in key:
            output_columns.append(key.replace("Allele 1: ", "Allele: "))
        elif "Allele 2: " in key:
            pass
        else:
            output_columns.append(key)

    output_rows = []
    for row in df.to_dict("records"):
        output_row1 = {}
        output_row2 = {}
        for key in df.columns:
            if "Allele 1: " in key:
                output_key = key.replace("Allele 1: ", "Allele: ")
                output_row1[output_key] = row[key]
            elif "Allele 2: " in key:
                output_key = key.replace("Allele 2: ", "Allele: ")
                output_row2[output_key] = row[key]
            else:
                output_key = key
                output_row1[output_key] = row[key]
                output_row2[output_key] = row[key]
        output_rows.append(output_row1)
        output_rows.append(output_row2)

    df = None

    alleles_df = pd.DataFrame(output_rows, columns=output_columns)
    alleles_df = alleles_df.astype(str).replace(r"\.0$", "", regex=True).replace(r"^nan$", "", regex=True).replace(r"^None$", "", regex=True)
    alleles_df.to_csv(output_tsv_path, index=False, header=True, sep="\t")

    print(f"Wrote {len(alleles_df):,d} rows to {output_tsv_path}")


if __name__ == "__main__":
    main()
