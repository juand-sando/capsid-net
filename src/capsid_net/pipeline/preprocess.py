import argparse
import os

import numpy as np
import pandas as pd


def rename_chains(df, rename_file):
    if rename_file and os.path.exists(rename_file):
        rename_df = pd.read_csv(rename_file)
        rename_map = dict(zip(rename_df["Chmx_rename"], rename_df["PISA_rename"]))
        df["Monomer1"] = df["Monomer1"].replace(rename_map)
        df["Monomer2"] = df["Monomer2"].replace(rename_map)
        print("Chain renaming applied.")
    else:
        print("No rename file provided, skipping renaming.")
    return df


def apply_grouping_from_file(df, grouping_file):
    if not os.path.exists(grouping_file):
        raise FileNotFoundError(f"Grouping file '{grouping_file}' not found.")

    grouping_df = pd.read_csv(grouping_file)
    protomer_to_group = dict(zip(grouping_df["Protomer"], grouping_df["Group"]))

    def map_to_group(protomer):
        if "_x" in protomer:
            base, suffix = protomer.split("_x", 1)
            group = protomer_to_group.get(base, base)
            return f"{group}_x{suffix}"
        return protomer_to_group.get(protomer, protomer)

    df["Monomer1"] = df["Monomer1"].apply(map_to_group)
    df["Monomer2"] = df["Monomer2"].apply(map_to_group)

    original_len = len(df)
    df = df[df["Monomer1"] != df["Monomer2"]].copy()
    removed = original_len - len(df)
    print(f"Grouping applied from {grouping_file}. Removed {removed} self-interactions.")
    return df


def add_redundant_interactions(df):
    group_cols = ["Monomer1", "Monomer2"]
    sum_cols = ["Area", "DeltaG", "Nhb", "Nsb", "Nds"]
    keep_cols = ["##", "Id", "Symmetry_operation", "Sym_Id"]

    df = df.sort_values(by=keep_cols)

    grouped_df = df.groupby(group_cols, as_index=False).agg(
        {**{col: "first" for col in keep_cols}, **{col: "sum" for col in sum_cols}}
    )

    print(f"Reduced redundant interactions: {len(df)} -> {len(grouped_df)} rows")
    return grouped_df


def create_interaction_matrix(df):
    proteins = sorted(set(df["Monomer1"]).union(set(df["Monomer2"])))
    protein_index = {protein: i for i, protein in enumerate(proteins)}

    interaction_matrix = np.zeros((len(proteins), len(proteins)))

    for _, row in df.iterrows():
        i = protein_index[row["Monomer1"]]
        j = protein_index[row["Monomer2"]]
        interaction_matrix[i, j] = row["Area"]
        interaction_matrix[j, i] = row["Area"]

    return interaction_matrix, proteins


def remove_self_interactions(df):
    original_len = len(df)
    filtered_df = df[df["Monomer1"] != df["Monomer2"]].copy()
    removed = original_len - len(filtered_df)
    if removed:
        print(f"Removed {removed} self-interactions before matrix construction.")
    return filtered_df


def run_analysis(args):
    os.makedirs(args.output, exist_ok=True)

    df = pd.read_csv(args.interactions)
    df = rename_chains(df, args.rename_file)

    if args.grouping_file:
        df = apply_grouping_from_file(df, args.grouping_file)
        df = add_redundant_interactions(df)
    else:
        print("No grouping file provided; skipping explicit regrouping.")

    df = remove_self_interactions(df)
    interaction_matrix, proteins = create_interaction_matrix(df)
    pd.DataFrame(df).to_csv(os.path.join(args.output, "pisa_processed.csv"))
    pd.DataFrame(interaction_matrix, index=proteins, columns=proteins).to_csv(os.path.join(args.output, "interaction_matrix.csv"))


def build_parser(parser):
    parser.add_argument("--interactions", "--i", dest="interactions", required=True, help="Path to the interactions CSV file")
    parser.add_argument("--output", required=True, help="Folder where results will be saved")
    parser.add_argument("--rename_file", default=None, help="Optional path to rename CSV")
    parser.add_argument("--grouping_file", help="Optional path to grouping CSV")
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Process interaction tables into grouped interaction matrices.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
