import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd

from capsid_net.utils.prepared_config import write_preprocess_config


def load_dataframes(interactions_file, com_file, prot_classes_file, rename_file, grouping_file):
    interactions_df = pd.read_csv(interactions_file)
    com_df = pd.read_csv(com_file)
    classes_df = pd.read_csv(prot_classes_file)
    rename_df = pd.read_csv(rename_file)
    grouping_df = pd.read_csv(grouping_file)
    return interactions_df, com_df, classes_df, rename_df, grouping_df


def build_group_map(grouping_df):
    return dict(zip(grouping_df["Protomer"], grouping_df["Group"]))


def map_group_label(label, protomer_to_group):
    if "_x" in label:
        base, suffix = label.split("_x", 1)
        group = protomer_to_group.get(base, base)
        return f"{group}_x{suffix}"
    return protomer_to_group.get(label, label)


def rename_interaction_chains(df, rename_df):
    rename_map = dict(zip(rename_df["Chmx_rename"], rename_df["PISA_rename"]))
    df["Monomer1"] = df["Monomer1"].replace(rename_map)
    df["Monomer2"] = df["Monomer2"].replace(rename_map)
    print("Interaction-chain renaming applied.")
    return df


def rename_com_chains(com_df, rename_df):
    rename_dict = dict(zip(rename_df["Chmx_rename"], rename_df["PISA_rename"]))
    com_df["chain"] = com_df["chain"].map(rename_dict).fillna(com_df["chain"])
    print("COM-chain renaming applied.")
    return com_df


def apply_grouping_to_interactions(df, grouping_df):
    protomer_to_group = build_group_map(grouping_df)

    df["Monomer1"] = df["Monomer1"].apply(lambda value: map_group_label(value, protomer_to_group))
    df["Monomer2"] = df["Monomer2"].apply(lambda value: map_group_label(value, protomer_to_group))

    original_len = len(df)
    df = df[df["Monomer1"] != df["Monomer2"]].copy()
    removed = original_len - len(df)
    print(f"Interaction grouping applied. Removed {removed} self-interactions.")
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


def remove_self_interactions(df):
    original_len = len(df)
    filtered_df = df[df["Monomer1"] != df["Monomer2"]].copy()
    removed = original_len - len(filtered_df)
    if removed:
        print(f"Removed {removed} self-interactions before matrix construction.")
    return filtered_df


def create_interaction_matrix(df):
    proteins = sorted(set(df["Monomer1"]).union(set(df["Monomer2"])))
    protein_index = {protein: i for i, protein in enumerate(proteins)}

    interaction_matrix = np.zeros((len(proteins), len(proteins)))

    for _, row in df.iterrows():
        i = protein_index[row["Monomer1"]]
        j = protein_index[row["Monomer2"]]
        interaction_matrix[i, j] = row["Area"]
        interaction_matrix[j, i] = row["Area"]

    return pd.DataFrame(interaction_matrix, index=proteins, columns=proteins)


def create_grouped_coms(com_df, grouping_df, classes_df):
    protomer_to_group = build_group_map(grouping_df)
    class_by_protein = dict(zip(classes_df["Protein_name"], classes_df["Class"]))

    grouped = {}
    for _, row in com_df.iterrows():
        chain_name = row["chain"]
        group_name = map_group_label(chain_name, protomer_to_group)
        base_name = chain_name.split("_", 1)[0]
        protein_class = str(class_by_protein.get(base_name, "")).strip()

        if group_name not in grouped:
            grouped[group_name] = {"coords": [], "classes": set()}

        grouped[group_name]["coords"].append(np.array([row["x"], row["y"], row["z"]], dtype=float))
        if protein_class:
            grouped[group_name]["classes"].add(protein_class)

    results = []
    for group_name, payload in grouped.items():
        center = np.stack(payload["coords"]).mean(axis=0)
        classes = sorted(payload["classes"])
        if len(classes) == 1:
            class_name = classes[0]
        elif len(classes) > 1:
            class_name = "Mixed"
        else:
            class_name = ""

        results.append(
            {
                "group": group_name,
                "rx": center[0],
                "ry": center[1],
                "rz": center[2],
                "Class": class_name,
            }
        )

    return pd.DataFrame(results).sort_values("group").reset_index(drop=True)


def filter_mcp_protomers(com_df, classes_df):
    mcp_basenames = classes_df[classes_df["Class"] == "MCP"]["Protein_name"].tolist()
    return com_df[com_df["chain"].apply(lambda chain: any(chain.startswith(name) for name in mcp_basenames))].copy()


def build_group_dictionaries(com_df, grouping_df):
    protomer_to_group = build_group_map(grouping_df)
    groupings = {}

    for _, row in com_df.iterrows():
        chain_name = row["chain"]
        base_name = chain_name.split("_", 1)[0]
        if base_name not in protomer_to_group:
            continue

        group_name = map_group_label(chain_name, protomer_to_group)
        groupings.setdefault(group_name, {})[chain_name] = np.array([row["x"], row["y"], row["z"]], dtype=float)

    return groupings


def get_center_and_normal(point_dict):
    points = list(point_dict.values())
    if len(points) != 3:
        raise ValueError("point_dict must contain exactly 3 points.")

    a, b, c = points
    center = (a + b + c) / 3
    u = b - a
    v = c - a
    normal = np.cross(u, v)

    if np.dot(normal, center) < 0:
        normal = -normal

    normal_unit = normal / np.linalg.norm(normal)
    return np.concatenate((center, normal_unit))


def compute_centers_and_normals(groupings):
    results = []
    for group_name, point_dict in groupings.items():
        if len(point_dict) == 3:
            result = get_center_and_normal(point_dict)
            results.append(
                {
                    "group": group_name,
                    "rx": result[0],
                    "ry": result[1],
                    "rz": result[2],
                    "nx": result[3],
                    "ny": result[4],
                    "nz": result[5],
                }
            )
        else:
            print(f"Skipping group '{group_name}' - it has {len(point_dict)} protomers (needs exactly 3).")
    return pd.DataFrame(results)


def create_capsid_angles(center_normal_df, interaction_df):
    group_names = center_normal_df["group"].values
    missing = [name for name in group_names if name not in interaction_df.index]
    if missing:
        raise ValueError(f"The following group names are missing from interaction_matrix.csv: {missing}")

    interaction_subset = interaction_df.loc[group_names, group_names].values
    interaction_mask = interaction_subset > 0

    normals = center_normal_df[["nx", "ny", "nz"]].values
    normals = normals / np.linalg.norm(normals, axis=1, keepdims=True)

    num = len(normals)
    angle_matrix = np.zeros((num, num))

    for i in range(num):
        for j in range(i + 1, num):
            if interaction_mask[i, j]:
                ni = normals[i]
                nj = normals[j]

                dot = np.clip(np.dot(ni, nj), -1.0, 1.0)
                angle = np.degrees(np.arccos(dot))

                if angle > 90:
                    dot_flipped = np.clip(np.dot(ni, -nj), -1.0, 1.0)
                    angle_flipped = np.degrees(np.arccos(dot_flipped))

                    if angle_flipped <= 90:
                        angle = angle_flipped
                    else:
                        raise ValueError(
                            f"Cannot resolve orientation between vectors {i} and {j} "
                            f"(angle={angle:.2f} deg, flipped={angle_flipped:.2f} deg)"
                        )

                angle_matrix[i, j] = angle
                angle_matrix[j, i] = angle

    return pd.DataFrame(angle_matrix, index=group_names, columns=group_names)


def run_analysis(args):
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    interactions_df, com_df, classes_df, rename_df, grouping_df = load_dataframes(
        args.interactions,
        args.com,
        args.prot_classes,
        args.rename,
        args.grouping,
    )

    interactions_df = rename_interaction_chains(interactions_df, rename_df)
    interactions_df = apply_grouping_to_interactions(interactions_df, grouping_df)
    interactions_df = add_redundant_interactions(interactions_df)
    interactions_df = remove_self_interactions(interactions_df)

    interaction_matrix_df = create_interaction_matrix(interactions_df)
    interactions_path = output_dir / "pisa_processed.csv"
    interaction_matrix_path = output_dir / "interaction_matrix.csv"
    interactions_df.to_csv(interactions_path, index=False)
    interaction_matrix_df.to_csv(interaction_matrix_path)

    com_df = rename_com_chains(com_df, rename_df)
    grouped_coms_df = create_grouped_coms(com_df, grouping_df, classes_df)
    grouped_coms_path = output_dir / "grouped_coms.csv"
    grouped_coms_df.to_csv(grouped_coms_path, index=False)

    mcp_com_df = filter_mcp_protomers(com_df, classes_df)
    mcp_groupings = build_group_dictionaries(mcp_com_df, grouping_df)
    center_normal_df = compute_centers_and_normals(mcp_groupings)
    capsid_angles_df = create_capsid_angles(center_normal_df, interaction_matrix_df)
    capsid_angles_path = output_dir / "capsid_angles.csv"
    capsid_angles_df.to_csv(capsid_angles_path)

    config_path = output_dir / "preprocess_config.json"
    write_preprocess_config(
        config_path,
        {
            "schema_version": 1,
            "output_dir": str(output_dir.resolve()),
            "inputs": {
                "interactions": str(Path(args.interactions).resolve()),
                "com": str(Path(args.com).resolve()),
                "prot_classes": str(Path(args.prot_classes).resolve()),
                "rename": str(Path(args.rename).resolve()),
                "grouping": str(Path(args.grouping).resolve()),
            },
            "outputs": {
                "pisa_processed": str(interactions_path.resolve()),
                "interaction_matrix": str(interaction_matrix_path.resolve()),
                "grouped_coms": str(grouped_coms_path.resolve()),
                "capsid_angles": str(capsid_angles_path.resolve()),
            },
        },
    )

    print(f"Saved processed interactions to: {interactions_path}")
    print(f"Saved interaction matrix to: {interaction_matrix_path}")
    print(f"Saved grouped COMs to: {grouped_coms_path}")
    print(f"Saved capsid angles to: {capsid_angles_path}")
    print(f"Saved preprocess config to: {config_path}")


def build_parser(parser):
    parser.add_argument("--interactions", "--i", dest="interactions", required=True, help="Path to the raw interactions CSV file")
    parser.add_argument("--com", required=True, help="Path to the chain center-of-mass CSV file")
    parser.add_argument("--prot_classes", required=True, help="Path to the protein classes CSV file")
    parser.add_argument("--rename", required=True, help="Path to the rename CSV file")
    parser.add_argument("--grouping", required=True, help="Path to the grouping CSV file")
    parser.add_argument("--output", required=True, help="Folder where prepared results will be saved")
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Preprocess a PISA interaction file and related metadata for downstream capsid plots.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
