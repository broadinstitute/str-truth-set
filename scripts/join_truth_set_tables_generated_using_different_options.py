"""This script left-joins two truth set tables generated using different parameters or command line args."""

import argparse
import pandas as pd

def parse_args_and_load_input_tables():
    parser = argparse.ArgumentParser()
    parser.add_argument("--table1", required=True)
    parser.add_argument("--table2", required=True)
    parser.add_argument("--table2-suffix", required=True, help="Suffix to append to column names in table2")
    parser.add_argument("--output-table", required=True)
    parser.add_argument("-k", "--key-column", action="append")
    parser.add_argument("-c", "--column-to-keep", action="append")
    parser.add_argument("--verbose", action="store_true")
    args = parser.parse_args()

    df1 = pd.read_table(args.table1)
    print(f"Parsed {len(df1):,d} rows from {args.table1}")
    df2 = pd.read_table(args.table2)
    print(f"Parsed {len(df2):,d} rows from {args.table2}")

    valid_columns = set(df1.columns) & set(df2.columns)

    if not args.key_column:
        parser.error("-k arg not specified.")
    invalid_key_columns = [k for k in args.key_column if k not in valid_columns]
    if invalid_key_columns:
        parser.error("Invalid column names: " + ", ".join(invalid_key_columns) + ". "
                     "Must be one or more columns from: " + ", ".join(valid_columns))

    if not args.column_to_keep:
        parser.error("-c arg not specified. Must be one or more columns from: " + ", ".join(valid_columns))
    invalid_columns_to_keep = [k for k in args.column_to_keep if k not in df1.columns or k not in df2.columns]
    if invalid_columns_to_keep:
        parser.error("Invalid column names: " + ", ".join(invalid_columns_to_keep) + ". "
                     "Must be one or more columns from: " + ", ".join(valid_columns))

    df2 = df2[list(args.key_column) + list(args.column_to_keep)]

    return args, df1, df2


def main():
    args, df1, df2 = parse_args_and_load_input_tables()

    df1 = df1.set_index(args.key_column)
    df2 = df2.set_index(args.key_column)

    print("Left join using key:", ", ".join(args.key_column))
    df = df1.join(df2, how="left", lsuffix="", rsuffix=args.table2_suffix).reset_index()

    for column_to_keep in args.column_to_keep:
        column_to_keep_with_suffix = column_to_keep + args.table2_suffix
        missing_values_count = sum(df[column_to_keep_with_suffix].isna())
        print(f"Filling in {missing_values_count:,d} out of {len(df)} ({100*missing_values_count/len(df):0.1}%) "
              f"missing values in {column_to_keep_with_suffix} from {column_to_keep}")
        df.loc[:, column_to_keep_with_suffix] = df[column_to_keep_with_suffix].fillna(df[column_to_keep])

    df.to_csv(args.output_table, sep="\t", index=False)

    if args.verbose:
        print("Output table has columns:\n", ", ".join(df.columns))
    print(f"Wrote {len(df):,d} rows to {args.output_table}")


if __name__ == "__main__":
    main()
