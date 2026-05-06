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

from capsid_net.utils.plot_styles import load_label_style_map

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


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


def plot_angle_network(
    center_normal_df,
    distance_matrix,
    angle_matrix,
    output_dir,
    title="Intercapsomer angle",
    tag_by_label=None,
    node_color_by_label=None,
    label_color_by_label=None,
):
    group_labels = center_normal_df["group"].values
    n = len(group_labels)

    points_3d = center_normal_df[["rx", "ry", "rz"]].values
    pca = PCA(n_components=2)
    points_2d = pca.fit_transform(points_3d)

    pos_dict = {group_labels[i]: points_2d[i] for i in range(len(group_labels))}

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

    node_colors = [node_color_by_label.get(label, "lightgrey") for label in group_labels]

    fig, ax = plt.subplots(figsize=(10, 10), constrained_layout=True)
    nx.draw_networkx_nodes(graph, pos_dict, node_color=node_colors, node_size=500, edgecolors="black", linewidths=0.8, ax=ax)
    label_dict = {label: tag_by_label.get(label, label) for label in group_labels}
    labels_by_color = {}
    for label in group_labels:
        color = label_color_by_label.get(label, "black")
        labels_by_color.setdefault(color, {})[label] = label_dict[label]

    for color, labels in labels_by_color.items():
        nx.draw_networkx_labels(
            graph,
            pos_dict,
            labels=labels,
            font_size=16,
            font_weight="bold",
            ax=ax,
            font_family="Arial",
            font_color=color,
        )
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
    tag_by_label = {}
    node_color_by_label = {}
    label_color_by_label = {}
    if args.tag_color_csv:
        tag_by_label, node_color_by_label, label_color_by_label = load_label_style_map(args.tag_color_csv)

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

    plot_angle_network(
        center_normal_df,
        distance_matrix,
        angle_matrix,
        args.output,
        tag_by_label=tag_by_label,
        node_color_by_label=node_color_by_label,
        label_color_by_label=label_color_by_label,
    )


def build_parser(parser):
    parser.add_argument("--com", required=True, help="Center of mass CSV file")
    parser.add_argument("--prot_classes", required=True, help="Protein class CSV file")
    parser.add_argument("--rename", required=True, help="CSV file with rename instructions")
    parser.add_argument("--grouping", required=True, help="CSV file with grouping information")
    parser.add_argument("--interaction", required=True, help="CSV file with interaction matrix")
    parser.add_argument("--custom_filter", help="Optional .txt file with group names to include in the network")
    parser.add_argument("--output", default="results", help="Output directory for CSVs and plots")
    parser.add_argument(
        "--tag_color_csv",
        help="Optional CSV with columns: label, final_tag, node_color, label_color",
    )
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Capsomer normal and distance analysis.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
