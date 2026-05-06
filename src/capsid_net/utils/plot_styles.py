import pandas as pd


def load_label_style_map(csv_path):
    mapping_df = pd.read_csv(csv_path)

    required_columns = {"label", "final_tag", "node_color", "label_color"}
    missing = required_columns - set(mapping_df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {csv_path}: {missing}")

    for column in ["label", "final_tag", "node_color", "label_color"]:
        mapping_df[column] = mapping_df[column].astype("string").str.strip()

    mapping_df = mapping_df[mapping_df["label"] != ""]

    tag_by_label = dict(zip(mapping_df["label"], mapping_df["final_tag"]))
    node_color_by_label = dict(zip(mapping_df["label"], mapping_df["node_color"]))
    label_color_by_label = dict(zip(mapping_df["label"], mapping_df["label_color"]))
    return tag_by_label, node_color_by_label, label_color_by_label
