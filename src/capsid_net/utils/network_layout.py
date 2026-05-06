import yaml
import numpy as np


def load_layout_config(path):
    if not path:
        return {"node_sizes": {}, "rules": []}

    with open(path, "r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}

    if not isinstance(data, dict):
        raise ValueError("Layout config must be a YAML mapping.")

    node_sizes = data.get("node_sizes", {}) or {}
    rules = data.get("rules", []) or []

    if not isinstance(node_sizes, dict):
        raise ValueError("'node_sizes' must be a mapping if provided.")
    if not isinstance(rules, list):
        raise ValueError("'rules' must be a list if provided.")

    return {"node_sizes": node_sizes, "rules": rules}


def compute_node_sizes(grouped_coms_df, layout_config):
    node_sizes_cfg = layout_config.get("node_sizes", {})
    default_size = node_sizes_cfg.get("default", 250)
    by_class = node_sizes_cfg.get("by_class", {}) or {}
    by_label = node_sizes_cfg.get("by_label", {}) or {}

    sizes = []
    for _, row in grouped_coms_df.iterrows():
        label = row["group"]
        class_name = row.get("Class", "")
        size = default_size

        if class_name in by_class:
            size = by_class[class_name]
        if label in by_label:
            size = by_label[label]

        sizes.append(size)

    return sizes


def apply_layout_rules(grouped_coms_df, points_2d, layout_config):
    rules = layout_config.get("rules", [])
    if not rules:
        return points_2d

    labels = grouped_coms_df["group"].tolist()
    classes = grouped_coms_df["Class"].tolist() if "Class" in grouped_coms_df.columns else [""] * len(grouped_coms_df)
    label_to_index = {label: idx for idx, label in enumerate(labels)}

    adjusted = points_2d.copy()

    for idx, rule in enumerate(rules, start=1):
        rule_type = rule.get("type")
        if rule_type == "label_shift":
            _apply_label_shift(adjusted, label_to_index, rule, idx)
        elif rule_type == "class_shift":
            _apply_class_shift(adjusted, classes, rule, idx)
        elif rule_type == "region_shift":
            _apply_region_shift(adjusted, label_to_index, rule, idx)
        else:
            raise ValueError(f"Unsupported network layout rule type: {rule_type!r}")

    return adjusted


def _apply_label_shift(points_2d, label_to_index, rule, rule_index):
    labels = rule.get("labels", [])
    dx = float(rule.get("dx", 0))
    dy = float(rule.get("dy", 0))

    for label in labels:
        if label not in label_to_index:
            print(f"Warning: layout rule {rule_index} references missing label '{label}'.")
            continue
        points_2d[label_to_index[label]] += np.array([dx, dy])


def _apply_class_shift(points_2d, classes, rule, rule_index):
    target_classes = set(rule.get("classes", []))
    dx = float(rule.get("dx", 0))
    dy = float(rule.get("dy", 0))

    if not target_classes:
        print(f"Warning: class_shift rule {rule_index} has no classes.")
        return

    for idx, class_name in enumerate(classes):
        if class_name in target_classes:
            points_2d[idx] += np.array([dx, dy])


def _apply_region_shift(points_2d, label_to_index, rule, rule_index):
    dx = float(rule.get("dx", 0))
    dy = float(rule.get("dy", 0))

    x_gt = _resolve_axis_threshold(points_2d, label_to_index, rule, axis=0, bound="gt", rule_index=rule_index)
    x_lt = _resolve_axis_threshold(points_2d, label_to_index, rule, axis=0, bound="lt", rule_index=rule_index)
    y_gt = _resolve_axis_threshold(points_2d, label_to_index, rule, axis=1, bound="gt", rule_index=rule_index)
    y_lt = _resolve_axis_threshold(points_2d, label_to_index, rule, axis=1, bound="lt", rule_index=rule_index)

    for idx, point in enumerate(points_2d):
        if x_gt is not None and not (point[0] > x_gt):
            continue
        if x_lt is not None and not (point[0] < x_lt):
            continue
        if y_gt is not None and not (point[1] > y_gt):
            continue
        if y_lt is not None and not (point[1] < y_lt):
            continue
        points_2d[idx] += np.array([dx, dy])


def _resolve_axis_threshold(points_2d, label_to_index, rule, axis, bound, rule_index):
    numeric_key = f"{'x' if axis == 0 else 'y'}_{bound}"
    label_key = f"{numeric_key}_label"

    if numeric_key in rule and label_key in rule:
        raise ValueError(f"Layout rule {rule_index} cannot define both '{numeric_key}' and '{label_key}'.")

    if numeric_key in rule:
        return float(rule[numeric_key])

    if label_key in rule:
        label = rule[label_key]
        if label not in label_to_index:
            print(f"Warning: layout rule {rule_index} references missing anchor label '{label}'.")
            return None
        return float(points_2d[label_to_index[label]][axis])

    return None
