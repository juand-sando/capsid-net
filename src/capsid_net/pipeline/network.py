import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import colormaps
from sklearn.decomposition import PCA

from capsid_net.utils.plot_styles import load_label_style_map

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def load_dataframes(interaction_file, com_file, prot_classes_file, rename_file, grouping_file):
    interaction_df = pd.read_csv(interaction_file, index_col=0)
    interaction_matrix = interaction_df.values
    interaction_labels = interaction_df.index.tolist()

    com_df = pd.read_csv(com_file)
    classes_df = pd.read_csv(prot_classes_file)
    rename_df = pd.read_csv(rename_file)
    grouping_df = pd.read_csv(grouping_file)
    return interaction_matrix, interaction_labels, com_df, classes_df, rename_df, grouping_df


def apply_renaming(com_df, rename_df):
    rename_dict = dict(zip(rename_df["Chmx_rename"], rename_df["PISA_rename"]))
    com_df["chain"] = com_df["chain"].map(rename_dict).fillna(com_df["chain"])
    return com_df


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


def get_centers(point_dict):
    points = list(point_dict.values())
    if len(points) < 1:
        raise ValueError("point_dict must contain at least one point.")
    points_array = np.stack(points)
    return points_array.mean(axis=0)


def compute_centers(groupings):
    results = []
    for group_name, point_dict in groupings.items():
        if len(point_dict) > 0:
            center = get_centers(point_dict)
            results.append({"group": group_name, "rx": center[0], "ry": center[1], "rz": center[2]})
        else:
            print(f"Skipping group '{group_name}' - it has no protomers.")
    return pd.DataFrame(results)


def eliminate_fastener_centers(center_df, classes_df):
    fasteners = set(classes_df[classes_df.iloc[:, 1] == "Fastener"].iloc[:, 0])
    hidden_fasteners = {f"{f}_x" for f in fasteners}

    def is_fastener_group(group):
        if group in fasteners:
            return True
        for hidden_fastener in hidden_fasteners:
            if group.startswith(hidden_fastener):
                return True
        return False

    filtered_df = center_df[~center_df["group"].apply(is_fastener_group)].copy()
    print(f"Removed {len(center_df) - len(filtered_df)} fastener centers.")
    return filtered_df


def summarize_centers(com_df, grouping_df, classes_df):
    groupings = build_group_dictionaries(com_df, grouping_df)
    center_df = compute_centers(groupings)

    grouped_chains = {chain for group in groupings.values() for chain in group}
    ungrouped_df = com_df[~com_df["chain"].isin(grouped_chains)].copy()
    ungrouped_df = ungrouped_df.rename(columns={"chain": "group", "x": "rx", "y": "ry", "z": "rz"})
    ungrouped_df = ungrouped_df[["group", "rx", "ry", "rz"]]

    all_centers_df = pd.concat([center_df, ungrouped_df], ignore_index=True)
    all_centers_df = eliminate_fastener_centers(all_centers_df, classes_df)
    return all_centers_df


def create_saltbridges_matrix(df):
    proteins = sorted(set(df["Monomer1"]).union(set(df["Monomer2"])))
    protein_index = {protein: i for i, protein in enumerate(proteins)}

    interaction_matrix = np.zeros((len(proteins), len(proteins)))

    for _, row in df.iterrows():
        i = protein_index[row["Monomer1"]]
        j = protein_index[row["Monomer2"]]
        interaction_matrix[i, j] = row["Dsb"]
        interaction_matrix[j, i] = row["Dsb"]

    return interaction_matrix, proteins


def plot_interaction_network(
    all_centers_df,
    interaction_matrix,
    interaction_labels,
    classes_df,
    output_dir,
    title="Surface Area Interactions",
    tag_by_label=None,
    node_color_by_label=None,
    label_color_by_label=None,
):
    all_centers_df = all_centers_df.set_index("group").reindex(interaction_labels).reset_index()

    valid_mask = ~all_centers_df["group"].astype("string").fillna("").str.startswith(("m", "n"))
    valid_labels = all_centers_df[valid_mask]["group"].tolist()

    label_to_index = {label: idx for idx, label in enumerate(interaction_labels)}
    valid_indices = [label_to_index[label] for label in valid_labels if label in label_to_index]

    interaction_matrix = interaction_matrix[np.ix_(valid_indices, valid_indices)]
    all_centers_df = all_centers_df.set_index("group").loc[valid_labels].reset_index()
    interaction_labels = valid_labels

    group_labels = all_centers_df["group"].values
    n = len(group_labels)
    coords = all_centers_df[["rx", "ry", "rz"]].to_numpy(dtype=float)

    pca = PCA(n_components=2, random_state=0)
    points_2d = pca.fit_transform(coords)

    group_list = list(all_centers_df["group"])

    if "pp_x2" in group_list and "tm" in group_list:
        pp_x2_index = group_list.index("pp_x2")
        tm_index = group_list.index("tm")
        tmx3_index = group_list.index("tm_x3")
        zazbx3_index = group_list.index("Za-Zb_x3")
        zczdx3_index = group_list.index("Zc-Zd_x3")
        zezfx3_index = group_list.index("Ze-Zf_x3")

        pp_x2_coords = points_2d[pp_x2_index]
        tm_coords = points_2d[tm_index]

        print(f"PCA coordinates for pp_x2: {pp_x2_coords}")
        print(f"PCA coordinates for tm: {tm_coords}")

        for i, point in enumerate(points_2d):
            if point[1] > tm_coords[1] and point[0] > pp_x2_coords[0]:
                points_2d[i] += np.array([100, 100])

        points_2d[tm_index] += np.array([30, 30])
        points_2d[tmx3_index] += np.array([-25, 0])
        points_2d[zazbx3_index] += np.array([-25, 0])
        points_2d[zczdx3_index] += np.array([-25, 0])
        points_2d[zezfx3_index] += np.array([-25, 0])
    else:
        print("Warning: 'pp_x2' or 'tm' not found in group labels.")

    delta_x, delta_y = 6.0, 7.0
    is_thread_mask = np.array([(label.startswith("r") and "_" not in label) for label in group_list], dtype=bool)
    points_2d[is_thread_mask] += np.array([delta_x, delta_y])
    is_thread_mask2 = np.array([(label.startswith("r") and "_x1" in label) for label in group_list], dtype=bool)
    points_2d[is_thread_mask2] += np.array([0, -delta_y])

    graph = nx.Graph()
    for label in group_labels:
        graph.add_node(label)

    edges = []
    edge_colors = []
    edge_widths = []

    cmap = colormaps["YlOrRd"]
    max_area = np.max(interaction_matrix)
    norm = plt.Normalize(vmin=0, vmax=max_area if max_area > 0 else 1)
    width_scale = 4

    for i in range(n):
        for j in range(i + 1, n):
            area = interaction_matrix[i, j]
            if area > 0:
                color = cmap(norm(area))
                edges.append((group_labels[i], group_labels[j]))
                edge_colors.append(color)
                edge_widths.append((area / max_area) * width_scale + 0.5)
                graph.add_edge(group_labels[i], group_labels[j])

    node_colors = [node_color_by_label.get(label, "#777777") for label in group_labels]

    fig, ax = plt.subplots(figsize=(10, 8), constrained_layout=True)

    class_map = dict(zip(classes_df["Protein_name"], classes_df["Class"]))
    node_sizes = []
    for label in group_labels:
        base_label = label.split("_x")[0]
        if class_map.get(base_label, "").lower() == "mcp":
            node_sizes.append(400)
        else:
            node_sizes.append(250)

    pos_dict = {group_labels[i]: points_2d[i] for i in range(n)}

    nx.draw_networkx_nodes(graph, pos_dict, node_color=node_colors, node_size=node_sizes, ax=ax, edgecolors="black", linewidths=0.8)

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
            font_family="Arial",
            font_size=16,
            font_color=color,
            font_weight="bold",
        )
    nx.draw_networkx_edges(graph, pos_dict, edgelist=edges, edge_color=edge_colors, width=edge_widths, ax=ax)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, orientation="horizontal")
    cbar.set_label("Interaction Buried Surface Area, BSA (A^2)", fontfamily="Arial", fontweight="bold", fontsize=13, labelpad=5)
    cbar.ax.tick_params(labelsize=12, width=1.5, length=10, pad=3, which="both", direction="out")

    ax.set_aspect("equal", adjustable="datalim")
    plt.axis("off")

    outputname = os.path.join(output_dir, "interaction_network_wo_SCPlabels_v2.svg")
    plt.savefig(outputname)
    plt.close()
    print(f"Network plot saved as {outputname}")


def run_analysis(args):
    os.makedirs(args.output, exist_ok=True)

    interaction_matrix, interaction_labels, com_df, classes_df, rename_df, grouping_df = load_dataframes(
        args.interaction, args.com, args.classes, args.rename, args.grouping
    )
    tag_by_label, node_color_by_label, label_color_by_label = load_label_style_map(args.tag_color_csv)

    if args.use_saltbridges:
        if not args.pisa_mod_file:
            raise ValueError("If --use_saltbridges is specified, --pisa_mod_file must be provided.")

        pisa_mod_df = pd.read_csv(args.pisa_mod_file)
        interaction_matrix, interaction_labels = create_saltbridges_matrix(pisa_mod_df)
        print("Loaded salt bridge interaction matrix (Dsb values)")
    else:
        print("Loaded area-based interaction matrix")

    com_df = apply_renaming(com_df, rename_df)
    all_centers_df = summarize_centers(com_df, grouping_df, classes_df)

    all_centers_path = os.path.join(args.output, "all_centers.csv")
    all_centers_df.to_csv(all_centers_path, index=False)
    print(f"All centers saved to {all_centers_path}")

    plot_interaction_network(
        all_centers_df,
        interaction_matrix,
        interaction_labels,
        classes_df,
        args.output,
        tag_by_label=tag_by_label,
        node_color_by_label=node_color_by_label,
        label_color_by_label=label_color_by_label,
    )


def build_parser(parser):
    parser.add_argument("--interaction", required=True, help="Path to interaction matrix CSV file")
    parser.add_argument("--com", required=True, help="Path to center of mass CSV file")
    parser.add_argument("--grouping", required=True, help="Path to groupings CSV file")
    parser.add_argument("--rename", required=True, help="Path to rename CSV file")
    parser.add_argument("--classes", required=True, help="Path to classes CSV file")
    parser.add_argument("--output", required=True, help="Output directory to save results")
    parser.add_argument("--use_saltbridges", action="store_true", help="Use Dsb values from a modified PISA CSV instead of area matrix")
    parser.add_argument("--pisa_mod_file", type=str, help="Path to PISA modified CSV containing Monomer1, Monomer2, and Dsb columns")
    parser.add_argument(
        "--tag_color_csv",
        type=str,
        required=True,
        help="CSV with columns: label, final_tag, node_color, label_color",
    )
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Analyze interaction network and generate a plot.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
