import argparse
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import LogNorm
from scipy.spatial.distance import cdist
from sklearn.decomposition import PCA

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def load_node_colors(node_colors_file):
    df = pd.read_csv(node_colors_file)
    cols = {c.lower(): c for c in df.columns}
    if "label" not in cols or "color" not in cols:
        raise ValueError("node_colors CSV must have columns 'label' and 'color'.")
    group_col = cols["label"]
    color_col = cols["color"]
    return dict(zip(df[group_col].astype(str), df[color_col].astype(str)))


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


def load_dataframes(com_file, prot_classes_file, rename_file, grouping_file):
    com_df = pd.read_csv(com_file)
    classes_df = pd.read_csv(prot_classes_file)
    rename_df = pd.read_csv(rename_file)
    grouping_df = pd.read_csv(grouping_file)
    return com_df, classes_df, rename_df, grouping_df


def apply_renaming(com_df, rename_df):
    rename_dict = dict(zip(rename_df["Chmx_rename"], rename_df["PISA_rename"]))
    com_df["chain"] = com_df["chain"].map(rename_dict).fillna(com_df["chain"])
    return com_df


def filter_mcp_protomers(com_df, classes_df):
    mcp_basenames = classes_df[classes_df["Class"] == "MCP"]["Protein_name"].tolist()
    return com_df[com_df["chain"].apply(lambda c: any(c.startswith(name) for name in mcp_basenames))].copy()


def build_group_dictionaries(capsomer_df, grouping_df):
    groupings = {}
    protomer_to_group = dict(zip(grouping_df["Protomer"], grouping_df["Group"]))

    for _, row in capsomer_df.iterrows():
        chain_name = row["chain"]
        base_name = chain_name.split("_")[0]
        suffix = chain_name[len(base_name):]

        if base_name in protomer_to_group:
            group_base = protomer_to_group[base_name]
            group_name = group_base + suffix
            if group_name not in groupings:
                groupings[group_name] = {}
            groupings[group_name][chain_name] = np.array([row["x"], row["y"], row["z"]])

    return groupings


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


def calculate_distances(df):
    coords = df[["rx", "ry", "rz"]].values
    return cdist(coords, coords)


def calculate_angles(df, interaction_file):
    interaction_df = pd.read_csv(interaction_file, index_col=0)
    group_names = df["group"].values

    missing = [name for name in group_names if name not in interaction_df.index]
    if missing:
        raise ValueError(f"The following group names are missing from interaction_matrix.csv: {missing}")

    interaction_subset = interaction_df.loc[group_names, group_names].values
    interaction_mask = interaction_subset > 0

    normals = df[["nx", "ny", "nz"]].values
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

    return angle_matrix


def plot_angle_network(center_normal_df, distance_matrix, angle_matrix, output_dir, title="Intercapsomer angle", node_colors_map=None):
    group_labels = center_normal_df["group"].values
    n = len(group_labels)

    points_3d = center_normal_df[["rx", "ry", "rz"]].values
    pca = PCA(n_components=2)
    points_2d = pca.fit_transform(points_3d)

    pos_dict = {group_labels[i]: points_2d[i] for i in range(len(group_labels))}

    had_subscript = {}
    display_label = {}
    for label in group_labels:
        if "_" in label:
            base = label.split("_", 1)[0]
            display_label[label] = base + "'"
            had_subscript[label] = True
        else:
            display_label[label] = label
            had_subscript[label] = False

    black_labels = {lab: display_label[lab] for lab in group_labels if "_" not in lab}
    white_labels = {lab: display_label[lab] for lab in group_labels if "_" in lab}

    graph = nx.Graph()
    for label in group_labels:
        graph.add_node(label)

    edges = []
    edge_colors = []
    edge_widths = []

    lower_limit = 1
    upper_limit = 30
    cmap = colormaps["YlOrRd"]
    norm = LogNorm(vmin=lower_limit, vmax=upper_limit)
    width = 5

    for i in range(n):
        for j in range(i + 1, n):
            angle = angle_matrix[i, j]
            if angle > 0:
                color = cmap(norm(angle))
                edges.append((group_labels[i], group_labels[j]))
                edge_colors.append(color)
                edge_widths.append(width)
                graph.add_edge(group_labels[i], group_labels[j])

    node_colors = []
    if node_colors_map:
        for label in group_labels:
            node_colors.append(node_colors_map.get(label, "lightgrey"))
    else:
        node_cmap = colormaps["bone"]
        min_frac = 0.40
        max_frac = 0.90
        subscript_numbers = [int(re.search(r"_x(\d+)", label).group(1)) for label in group_labels if "_x" in label]
        max_num = max(subscript_numbers) if subscript_numbers else 1
        for label in group_labels:
            match = re.search(r"_x(\d+)", label)
            if match:
                number = int(match.group(1))
                norm_val = (number / max_num) * (max_frac - min_frac) + min_frac
                node_colors.append(node_cmap(norm_val))
            else:
                node_colors.append("lightgrey")

    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    nx.draw_networkx_nodes(graph, pos_dict, node_color=node_colors, node_size=500, edgecolors="black", linewidths=0.8, ax=ax)
    for i, label in enumerate(group_labels):
        if had_subscript[label]:
            node_colors[i] = "white"

    nx.draw_networkx_labels(graph, pos_dict, labels=black_labels, font_size=16, font_weight="bold", ax=ax, font_family="Arial")
    nx.draw_networkx_labels(graph, pos_dict, labels=white_labels, font_size=16, font_weight="bold", ax=ax, font_family="Arial", font_color="white")
    nx.draw_networkx_edges(graph, pos_dict, edgelist=edges, edge_color=edge_colors, width=edge_widths, ax=ax)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, orientation="horizontal", extend="both")
    cbar.set_label("Inter-capsomer angle, phi (deg)", fontsize=13, fontweight="bold", labelpad=5)
    cbar.ax.tick_params(labelsize=12, width=1.5, length=10, pad=3, which="both", direction="out")

    plt.title(title)
    ax.set_aspect("equal", adjustable="datalim")
    plt.axis("off")

    plt.savefig(os.path.join(output_dir, "capsomer_angles.svg"), transparent=True)
    plt.close()
    print("Network plot saved as 'capsomer_angles.svg'")


def run_analysis(args):
    os.makedirs(args.output, exist_ok=True)

    com_df, classes_df, rename_df, grouping_df = load_dataframes(args.com, args.prot_classes, args.rename, args.grouping)
    node_colors_map = load_node_colors(args.node_colors) if args.node_colors else None

    com_df = apply_renaming(com_df, rename_df)
    capsomer_df = filter_mcp_protomers(com_df, classes_df)
    groupings_dict = build_group_dictionaries(capsomer_df, grouping_df)
    center_normal_df = compute_centers_and_normals(groupings_dict)

    if args.custom_filter:
        with open(args.custom_filter, "r") as handle:
            allowed_groups = set(line.strip() for line in handle if line.strip())
        before = len(center_normal_df)
        center_normal_df = center_normal_df[center_normal_df["group"].isin(allowed_groups)].reset_index(drop=True)
        after = len(center_normal_df)
        print(f"Applied custom filter: {after}/{before} groups retained")

    output_csv = os.path.join(args.output, "centers_and_normals.csv")
    center_normal_df.to_csv(output_csv, index=False)
    print(f"Saved center + normal vectors to: {output_csv}")

    distance_matrix = calculate_distances(center_normal_df)
    angle_matrix = calculate_angles(center_normal_df, args.interaction)

    dist_path = os.path.join(args.output, "capsomer_distances.csv")
    angle_path = os.path.join(args.output, "capsomer_angles.csv")

    pd.DataFrame(distance_matrix, index=center_normal_df["group"], columns=center_normal_df["group"]).to_csv(dist_path)
    pd.DataFrame(angle_matrix, index=center_normal_df["group"], columns=center_normal_df["group"]).to_csv(angle_path)

    print(f"Distance matrix saved to: {dist_path}")
    print(f"Angle matrix saved to: {angle_path}")

    plot_angle_network(center_normal_df, distance_matrix, angle_matrix, args.output, node_colors_map=node_colors_map)


def build_parser(parser):
    parser.add_argument("--com", required=True, help="Center of mass CSV file")
    parser.add_argument("--prot_classes", required=True, help="Protein class CSV file")
    parser.add_argument("--rename", required=True, help="CSV file with rename instructions")
    parser.add_argument("--grouping", required=True, help="CSV file with grouping information")
    parser.add_argument("--interaction", required=True, help="CSV file with interaction matrix")
    parser.add_argument("--custom_filter", help="Optional .txt file with group names to include in the network")
    parser.add_argument("--output", default="results", help="Output directory for CSVs and plots")
    parser.add_argument("--node_colors", help="Optional CSV mapping group label to color")
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Capsomer normal and distance analysis.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
