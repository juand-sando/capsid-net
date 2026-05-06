# capsid-net

CLI for preprocessing a PISA interaction file and analyzing capsid interaction networks.

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

Legacy wrapper scripts are also available:

```bash
python3 get_com.py --help
python3 preprocess.py --help
python3 angles.py --help
python3 network_builder_v3.py --help
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

Required columns:

```csv
Chmx_rename,PISA_rename
```

### `groupings.csv`

Required columns:

```csv
Protomer,Group
```

### `prot_classes.csv`

Required columns:

```csv
Protein_name,Class
```

### `tag_color.csv`

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

The exclusion file is plain text and supports two optional sections:

```text
Chains:
ABC_x2
TM_x1

Class:
Fastener
Decoration
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
