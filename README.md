# capsid-net

Initial CLI packaging refactor for capsid structure and PISA interaction analysis.

Install locally:

```bash
pip install -e .
```

The `com` stage now prefers a Gemmi-based parser for PDB and mmCIF input, with a legacy manual parser still available for comparison.

Available commands:

```bash
capsid-net --help
capsid-net com --help
capsid-net preprocess --help
capsid-net normals --help
capsid-net network --help
capsid-net run --help
```
