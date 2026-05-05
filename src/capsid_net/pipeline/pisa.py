import argparse
import os

import matplotlib
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
import seaborn as sns
from ete3 import NodeStyle
from ete3 import TextFace
from ete3 import Tree
from ete3 import TreeStyle
from scipy.spatial.distance import squareform


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


def group_similars(df, group_id, label):
    group_df = df[df["Id"] == group_id]
    unique_chains = set(group_df["Monomer1"]).union(set(group_df["Monomer2"]))
    group_map = {}

    for chain in unique_chains:
        interacting = set(group_df[group_df["Monomer1"] == chain]["Monomer2"]).union(
            group_df[group_df["Monomer2"] == chain]["Monomer1"]
        )
        interacting.add(chain)

        if any(c in group_map for c in interacting):
            continue

        base_names = {c[0] for c in interacting}
        suffixes = {c.split("_x")[-1] if "_x" in c else "" for c in interacting}

        if len(base_names) > 1:
            raise ValueError(f"[{label}] Conflictimg base names in group: {base_names}")
        if len(suffixes) > 1:
            raise ValueError(f"[{label}] Conflicting suffixes in the group: {base_names}")

        if group_id == 1:
            prefix = list(base_names)[0]
        if group_id == 2:
            prefix_list = [sorted(interacting)[0].split("_x")[0]]

        rep = f"{prefix}{'_x' + list(suffixes)[0] if list(suffixes)[0] else ''}"
        for c in interacting:
            group_map[c] = rep

    df.loc[df["Id"] != group_id, "Monomer1"] = df.loc[df["Id"] != group_id, "Monomer1"].replace(group_map)
    df.loc[df["Id"] != group_id, "Monomer2"] = df.loc[df["Id"] != group_id, "Monomer2"].replace(group_map)

    print(f"[{label}] Grouped {len(group_map)} chains.")
    return df[df["Id"] != group_id]


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


def eliminate_fasteners(df, class_csv):
    protein_classes = pd.read_csv(class_csv)
    fasteners = set(protein_classes[protein_classes.iloc[:, 1] == "Fastener"].iloc[:, 0])
    hidden_fasteners = {f"{f}_x" for f in fasteners}

    def is_fastener(monomer):
        if monomer in fasteners:
            return True
        for hidden_fastener in hidden_fasteners:
            if monomer.startswith(hidden_fastener):
                return True
        return False

    filtered_df = df[~df["Monomer1"].apply(is_fastener) & ~df["Monomer2"].apply(is_fastener)]
    print(f"Eliminated {len(df) - len(filtered_df)} fastener interactions.")
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

    return interaction_matrix, proteins


def create_relatedness_matrix(interaction_matrix):
    max_area = np.max(interaction_matrix)
    relatedness_matrix = 1 - (interaction_matrix / (max_area + 1))
    np.fill_diagonal(relatedness_matrix, 0)
    return relatedness_matrix


def generate_linkage_matrix(relatedness_matrix):
    return sch.linkage(squareform(relatedness_matrix), method="average")


def load_protein_classes(class_csv):
    df = pd.read_csv(class_csv)
    return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))


def assign_colors_to_classes(class_dict):
    unique_classes = set(class_dict.values())
    unique_classes.add("Neighbor")
    color_palette = sns.color_palette("husl", len(unique_classes))
    return {cls: matplotlib.colors.to_hex(color) for cls, color in zip(unique_classes, color_palette)}


def build_ete_tree(linkage_matrix, labels):
    n = len(labels)
    tree = {i: Tree(f"({labels[i]});") for i in range(n)}

    for i, (c1, c2, _, _) in enumerate(linkage_matrix, start=n):
        node = Tree()
        node.add_child(tree[int(c1)])
        node.add_child(tree[int(c2)])
        tree[i] = node

    return tree[max(tree.keys())]


def plot_circular_dendrogram(linkage_matrix, labels, class_csv, output_file):
    protein_classes = load_protein_classes(class_csv)
    class_colors = assign_colors_to_classes(protein_classes)
    tree = build_ete_tree(linkage_matrix, labels)

    tree_style = TreeStyle()
    tree_style.mode = "c"
    tree_style.show_leaf_name = False
    tree_style.scale = 100
    tree_style.draw_guiding_lines = False

    for node in tree.traverse():
        if node.is_leaf():
            protein_name = node.name
            protein_class = protein_classes.get(protein_name, "Neighbor")
            node_color = class_colors[protein_class]

            style = NodeStyle()
            style["size"] = 10
            style["fgcolor"] = node_color
            node.set_style(style)
            node.add_face(TextFace(protein_name, fsize=16, fgcolor=node_color), column=0, position="branch-right")

    tree.render(output_file, w=3000, dpi=300, tree_style=tree_style)
    print(f"Saved high-resolution circular dendrogram as {output_file}")


def recursive_split(linkage_matrix, labels, output_file, cluster_name="cluster"):
    num_items = len(labels)

    with open(output_file, "w") as handle:
        handle.write("Hierarchical Clustering Breakdown\n\n")

        def split_cluster(node_index, branch_name):
            if node_index < num_items:
                return [labels[node_index]]

            c1 = int(linkage_matrix[node_index - num_items, 0])
            c2 = int(linkage_matrix[node_index - num_items, 1])

            cluster1 = split_cluster(c1, f"{branch_name}.1")
            cluster2 = split_cluster(c2, f"{branch_name}.2")

            handle.write(f"{branch_name}.1:\n")
            handle.write(",".join(cluster1) + "\n")
            handle.write(f"{branch_name}.2:\n")
            handle.write(",".join(cluster2) + "\n\n")

            return cluster1 + cluster2

        root_index = len(linkage_matrix) + num_items - 1
        split_cluster(root_index, cluster_name)


def run_analysis(args):
    os.makedirs(args.output, exist_ok=True)

    df = pd.read_csv(args.interactions)
    df = rename_chains(df, args.rename_file)

    if args.grouping_file:
        df = apply_grouping_from_file(df, args.grouping_file)
        df = eliminate_fasteners(df, args.classes)
        df = add_redundant_interactions(df)

    grouped = False
    grouped_types = []

    if args.group_capsomers:
        print("Grouping capsomers...")
        df = group_similars(df, group_id=1, label="Capsomers")
        df = eliminate_fasteners(df, args.classes)
        df = add_redundant_interactions(df)
        grouped = True
        grouped_types.append("capsomers")

    if args.group_zippers:
        print("Grouping zippers...")
        df = group_similars(df, group_id=2, label="Zippers")
        df = add_redundant_interactions(df)
        grouped = True
        grouped_types.append("zippers")

    if grouped:
        if len(grouped_types) == 2:
            filename = "grouped_all.csv"
        elif grouped_types[0] == "capsomers":
            filename = "grouped_capsomers.csv"
        else:
            filename = "grouped_zippers.csv"
        grouped_path = os.path.join(args.output, filename)
        df.to_csv(grouped_path, index=False)
        print(f"Saved grouped dataframe as {grouped_path}")
    else:
        print("Skipped grouping!")

    interaction_matrix, proteins = create_interaction_matrix(df)
    relatedness_matrix = create_relatedness_matrix(interaction_matrix)
    linkage_matrix = generate_linkage_matrix(relatedness_matrix)

    pd.DataFrame(df).to_csv(os.path.join(args.output, "pisa_manipulated.csv"))
    pd.DataFrame(interaction_matrix, index=proteins, columns=proteins).to_csv(os.path.join(args.output, "interaction_matrix.csv"))
    pd.DataFrame(relatedness_matrix, index=proteins, columns=proteins).to_csv(os.path.join(args.output, "relatedness_matrix.csv"))

    recursive_split(linkage_matrix, proteins, os.path.join(args.output, "clusters.txt"))
    plot_circular_dendrogram(linkage_matrix, proteins, args.classes, os.path.join(args.output, "dendrogram.pdf"))


def build_parser(parser):
    parser.add_argument("--interactions", "--i", dest="interactions", required=True, help="Path to the interactions.csv file")
    parser.add_argument("--classes", required=True, help="Path to the protein_classes.csv file")
    parser.add_argument("--output", required=True, help="Folder where results will be saved")
    parser.add_argument("--rename_file", default=None, help="Optional path to rename CSV")
    parser.add_argument("--grouping_file", help="Optional path to grouping CSV")
    parser.add_argument("--group_capsomers", action="store_true", help="Enable capsomer grouping based on id=1 interactions")
    parser.add_argument("--group_zippers", action="store_true", help="Enable zipper grouping based on id=2 interactions")
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Generate hierarchical clusters and circular dendrogram.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
