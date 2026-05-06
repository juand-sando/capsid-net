import argparse
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import colormaps
from sklearn.decomposition import PCA

from capsid_net.utils.plot_filters import apply_plot_exclusions, load_plot_exclusions
from capsid_net.utils.plot_styles import load_label_style_map
from capsid_net.utils.network_layout import apply_layout_rules, compute_node_sizes, load_layout_config
from capsid_net.utils.prepared_config import load_preprocess_config

mpl.rcParams["font.family"] = "Arial"
mpl.rcParams["svg.fonttype"] = "none"


def load_prepared_inputs(config_file):
    config = load_preprocess_config(config_file)
    grouped_coms_df = pd.read_csv(config["outputs"]["grouped_coms"])
    interaction_df = pd.read_csv(config["outputs"]["interaction_matrix"], index_col=0)
    return config, grouped_coms_df, interaction_df


def filter_network_inputs(grouped_coms_df, interaction_df, exclude_file=None):
    grouped_coms_df = grouped_coms_df[grouped_coms_df["group"].isin(interaction_df.index)].copy()

    if exclude_file:
        exclusions = load_plot_exclusions(exclude_file)
        before = len(grouped_coms_df)
        grouped_coms_df = apply_plot_exclusions(grouped_coms_df, exclusions)
        after = len(grouped_coms_df)
        print(f"Applied plot exclusions: {after}/{before} groups retained")

    keep_labels = grouped_coms_df["group"].tolist()
    interaction_df = interaction_df.loc[keep_labels, keep_labels]
    return grouped_coms_df.reset_index(drop=True), interaction_df


def warn_on_missing_style_labels(group_labels, tag_by_label, node_color_by_label, label_color_by_label):
    missing = []
    for label in group_labels:
        if label not in tag_by_label or label not in node_color_by_label or label not in label_color_by_label:
            missing.append(label)

    if missing:
        print("Warning: the following plotted labels are missing one or more style entries in tag_color_csv:")
        for label in missing:
            print(f"  - {label}")


def plot_interaction_network(
    grouped_coms_df,
    interaction_df,
    output_dir,
    graph_format="svg",
    layout_config=None,
    title="Surface Area Interactions",
    tag_by_label=None,
    node_color_by_label=None,
    label_color_by_label=None,
):
    if len(grouped_coms_df) == 0:
        raise ValueError("No groups available to plot after applying filters.")

    group_labels = grouped_coms_df["group"].values
    warn_on_missing_style_labels(group_labels, tag_by_label, node_color_by_label, label_color_by_label)
    coords = grouped_coms_df[["rx", "ry", "rz"]].to_numpy(dtype=float)
    points_2d = PCA(n_components=2, random_state=0).fit_transform(coords)
    points_2d = apply_layout_rules(grouped_coms_df, points_2d, layout_config or {"node_sizes": {}, "rules": []})

    pos_dict = {group_labels[i]: points_2d[i] for i in range(len(group_labels))}
    graph = nx.Graph()
    for label in group_labels:
        graph.add_node(label)

    edges = []
    edge_colors = []
    edge_widths = []

    interaction_matrix = interaction_df.values
    cmap = colormaps["YlOrRd"]
    max_area = np.max(interaction_matrix)
    norm = plt.Normalize(vmin=0, vmax=max_area if max_area > 0 else 1)
    width_scale = 4

    n = len(group_labels)
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
    node_sizes = compute_node_sizes(grouped_coms_df, layout_config or {"node_sizes": {}, "rules": []})
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

    output_path = os.path.join(output_dir, f"capsid_interaction_network.{graph_format}")
    plt.savefig(output_path)
    plt.close()
    print(f"Network plot saved as {output_path}")


def run_analysis(args):
    config, grouped_coms_df, interaction_df = load_prepared_inputs(args.config)
    output_dir = args.output or config["output_dir"]
    os.makedirs(output_dir, exist_ok=True)

    grouped_coms_df, interaction_df = filter_network_inputs(grouped_coms_df, interaction_df, exclude_file=args.exclude)
    tag_by_label, node_color_by_label, label_color_by_label = load_label_style_map(args.tag_color_csv)
    layout_config = load_layout_config(args.layout_config)

    plot_interaction_network(
        grouped_coms_df,
        interaction_df,
        output_dir,
        graph_format=args.graph_format,
        layout_config=layout_config,
        tag_by_label=tag_by_label,
        node_color_by_label=node_color_by_label,
        label_color_by_label=label_color_by_label,
    )


def build_parser(parser):
    parser.add_argument("--config", required=True, help="Path to preprocess_config.json")
    parser.add_argument("--tag_color_csv", type=str, required=True, help="CSV with columns: label, final_tag, node_color, label_color")
    parser.add_argument("--exclude", help="Optional text file with exclusion sections 'Chains' and/or 'Class'")
    parser.add_argument("--layout-config", help="Optional YAML file with node size and layout adjustment rules")
    parser.add_argument("--output", help="Optional output directory override. Defaults to the preprocess output directory.")
    parser.add_argument("--graph-format", choices=["svg", "png"], default="svg", help="Output format for the network figure")
    parser.set_defaults(func=run_analysis)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot capsid interaction networks from preprocessed capsid data.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
