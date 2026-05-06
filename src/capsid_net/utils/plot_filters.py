def load_plot_exclusions(path):
    rules = {"Chains": set(), "Class": set()}
    current_section = None

    with open(path, "r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            header = line[:-1] if line.endswith(":") else line
            if header in rules:
                current_section = header
                continue

            if current_section is None:
                raise ValueError("Exclusion entries must be placed under a 'Chains' or 'Class' section.")

            rules[current_section].add(line)

    return rules


def apply_plot_exclusions(grouped_coms_df, exclusions):
    filtered_df = grouped_coms_df.copy()

    if exclusions["Chains"]:
        filtered_df = filtered_df[~filtered_df["group"].isin(exclusions["Chains"])]

    if exclusions["Class"] and "Class" in filtered_df.columns:
        filtered_df = filtered_df[~filtered_df["Class"].isin(exclusions["Class"])]

    return filtered_df.reset_index(drop=True)
