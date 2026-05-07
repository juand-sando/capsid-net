# capsid-net

CLI for processing a capsid model (PDB/CIF) and its associated PISA interaction file to analyze capsids through interaction networks.

## Overview

`capsid-net` builds interaction-network representations for complex viral capsids, including the large capsids found in Nucleocytoviricota.

The program is designed to take:

- a single structure file in `.pdb` or `.cif` format
  - this structure can contain just the asymmetric unit (ASU)
  - or the ASU plus surrounding copies if you want to inspect interactions at the ASU border
- the output of Protein Interfaces, Surfaces and Assemblies (PISA) in a clean CSV format (<https://www.ebi.ac.uk/pdbe/pisa/>)
  - the PISA table should correspond to the same structure that was used as the coordinate input

It then processes those inputs into reusable intermediate files and final network-style visualizations for downstream inspection.

The current workflow supports two closely related analysis modes:

### 1. Interaction-network analysis

`capsid-net` uses pairwise interface tables derived from PISA to build an outer capsid interaction network. In this workflow, buried surface area (BSA) values from individual protomer–protomer contacts are parsed, renamed, regrouped into higher-order functional units, and then aggregated into unit-level interaction strengths. This is intended to simplify the representation of large capsids and better reflect assembly-level organization rather than treating every chain-chain contact as an independent edge.

In practical terms, the preprocessing stage:

- reads the raw PISA interaction table
- applies the renaming and grouping scheme you provide
- sums interface measurements over all protomer–protomer contacts connecting the same pair of functional units
- writes a processed interaction table and an interaction matrix that can be rendered as a capsid interaction network

This is the part of the program used to study BSA-driven interaction patterns between the components of the capsid.

### 2. Capsid curvature analysis

`capsid-net` also supports curvature-style analysis focused on major capsid protein (MCP) capsomers. For this mode, the tool starts from the atomic coordinates in the input structure, computes the centre of mass of each protomer, and then uses the three MCP protomer centres associated with each capsomer to define a local surface patch and its outward normal vector.

Neighbouring MCP capsomers are then compared by measuring the angle between their surface normals. These pairwise capsomer angles are written to a matrix and can be rendered as an angle-colored graph for visual inspection of local curvature relationships across the capsid surface.

This is the part of the program used to study intercapsomer curvature geometry rather than interface strength alone.

## Install

### Conda environment

From the repository root, the simplest setup is:

```bash
conda env create -f environment.yml
conda activate capsid-net
```

This uses [`environment.yml`](/Users/juand.sando/Documents/fv3_repositories/capsid-net/environment.yml:1), which creates a Conda environment with Python and `pip`, then installs `capsid-net` in editable mode from the local checkout. The package dependencies are taken from [`pyproject.toml`](/Users/juand.sando/Documents/fv3_repositories/capsid-net/pyproject.toml:1).

After activation, verify the CLI is available:

```bash
capsid-net --help
```

### Manual install

```bash
conda create -n capsid-net python=3.11 pip -y
conda activate capsid-net
pip install -e .
```

## Commands

```bash
capsid-net --help
capsid-net com --help
capsid-net preprocess --help
capsid-net angles --help
capsid-net network --help
capsid-net run --help
```

## Workflow

The current workflow is:

1. `com`
   - read a structure file (`.pdb` or `.cif`)
   - write `center_of_mass.csv`
2. `preprocess`
   - read the raw PISA interaction CSV plus COM/class/rename/grouping metadata
   - write prepared files for both plots
3. `angles`
   - read the preprocess manifest and render the angle plot
4. `network`
   - read the preprocess manifest and render the interaction network plot

### End-to-end

```bash
capsid-net run \
  --structure my_capsid.cif \
  --interactions my_interactions.csv \
  --classes prot_classes.csv \
  --rename rename.csv \
  --grouping groupings.csv \
  --tag_color_csv tag_color.csv \
  --output results \
  --angles-exclude angles_exclude.txt \
  --network-exclude network_exclude.txt \
  --network-layout-config network_layout.yaml \
  --graph-format svg
```

### Step-by-step

```bash
capsid-net com \
  --structure my_capsid.cif \
  --output results/center_of_mass.csv
```

```bash
capsid-net preprocess \
  --interactions my_interactions.csv \
  --com results/center_of_mass.csv \
  --prot_classes prot_classes.csv \
  --rename rename.csv \
  --grouping groupings.csv \
  --output results
```

```bash
capsid-net angles \
  --config results/preprocess_config.json \
  --tag_color_csv tag_color.csv \
  --exclude angles_exclude.txt \
  --graph-format svg
```

```bash
capsid-net network \
  --config results/preprocess_config.json \
  --tag_color_csv tag_color.csv \
  --exclude network_exclude.txt \
  --layout-config network_layout.yaml \
  --graph-format svg
```

## Preprocess outputs

`preprocess` currently writes:

- `pisa_processed.csv`
- `interaction_matrix.csv`
- `grouped_coms.csv`
- `capsid_angles.csv`
- `preprocess_config.json`

`preprocess_config.json` is the handoff file used by `angles` and `network`. It contains:

- original input paths
- prepared output paths
- the default output directory

Example shape:

```json
{
  "schema_version": 1,
  "output_dir": "/abs/path/to/results",
  "inputs": {
    "interactions": "/abs/path/to/interactions.csv",
    "com": "/abs/path/to/center_of_mass.csv",
    "prot_classes": "/abs/path/to/prot_classes.csv",
    "rename": "/abs/path/to/rename.csv",
    "grouping": "/abs/path/to/groupings.csv"
  },
  "outputs": {
    "pisa_processed": "/abs/path/to/results/pisa_processed.csv",
    "interaction_matrix": "/abs/path/to/results/interaction_matrix.csv",
    "grouped_coms": "/abs/path/to/results/grouped_coms.csv",
    "capsid_angles": "/abs/path/to/results/capsid_angles.csv"
  }
}
```

## File contracts

### `rename.csv`

This file handles the fact that, when studying the interactions of an ASU plus its surroundings, equivalent chains in different ASUs can have different chain IDs. Here you can rename chain IDs to a shared base ID plus a suffix indicating which ASU they come from. This allows you to apply the same downstream rules across all copies. Renaming is also the first transformation applied, so subsequent files should use the renamed chain IDs.

Required columns:

```csv
Chmx_rename,PISA_rename
```

### `groupings.csv`

If you know that different chains belong to a specific unit, for example because the BSA of their interface is much higher than the rest, as in MCP capsomers, you can group them together. If you have renamed the chains of neighbouring ASUs to match the IDs of the original ASU plus their suffixes, you can apply one grouping logic to all of them. Otherwise, you may need to specify the logic for each group separately.

Required columns:

```csv
Protomer,Group
```

### `prot_classes.csv`

For each chain or grouped unit, define what entity it is, for example `MCP`, `TmP`, or `Penton`, so this annotation can be used downstream. Do not use whitespace in class names.

Required columns:

```csv
Protein_name,Class
```

### `tag_color.csv`

This file defines the displayed label and colors for nodes in the plots. Each label can be configured individually.

Required columns:

```csv
label,final_tag,node_color,label_color
```

This file controls:

- displayed node label text
- node fill color
- label text color

If a plotted label is missing one or more style entries, the program will print a warning and fall back to default styling.

### Exclusion file

The exclusion file is a plain-text list of nodes or classes to ignore during plotting, and it supports two optional sections. For example:

```text
Chains:
tm

Class:
TmP
```

Rules:

- entries under `Chains` match the `group` column in `grouped_coms.csv`
- entries under `Class` match the `Class` column in `grouped_coms.csv`
- blank lines are allowed
- lines starting with `#` are ignored

### Network layout YAML

This file is optional and only used by `network`.

Supported top-level keys:

- `node_sizes`
- `rules`

Example:

```yaml
node_sizes:
  default: 250
  by_class:
    MCP: 400
  by_label:
    tm: 320

rules:
  - type: label_shift
    labels: [tm]
    dx: 30
    dy: 30

  - type: class_shift
    classes: [Thread]
    dx: 6
    dy: 7

  - type: region_shift
    x_gt_label: pp_x2
    y_gt_label: tm
    dx: 100
    dy: 100
```

Supported rule types:

- `label_shift`
- `class_shift`
- `region_shift`

For `region_shift`, you can use either anchor labels or numeric thresholds:

- `x_gt_label`, `x_lt_label`, `y_gt_label`, `y_lt_label`
- `x_gt`, `x_lt`, `y_gt`, `y_lt`

Rules are applied in order.

## Examples directory

The repository includes toy example files in [`examples/`](./examples):

- [`examples/interactions.csv`](./examples/interactions.csv)
- [`examples/center_of_mass.csv`](./examples/center_of_mass.csv)
- [`examples/prot_classes.csv`](./examples/prot_classes.csv)
- [`examples/rename.csv`](./examples/rename.csv)
- [`examples/groupings.csv`](./examples/groupings.csv)
- [`examples/tag_color.csv`](./examples/tag_color.csv)
- [`examples/angles_exclude.txt`](./examples/angles_exclude.txt)
- [`examples/network_exclude.txt`](./examples/network_exclude.txt)
- [`examples/network_layout.yaml`](./examples/network_layout.yaml)
- [`examples/preprocess_config.example.json`](./examples/preprocess_config.example.json)

These are meant to illustrate file formats and CLI wiring. They are not intended as a scientifically meaningful reference dataset.
