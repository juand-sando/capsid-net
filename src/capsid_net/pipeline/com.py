import argparse
import csv
from pathlib import Path

try:
    import gemmi
except ImportError:  # pragma: no cover - depends on local environment
    gemmi = None


ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
}


def calculate_center_of_mass_legacy(structure_filename):
    chain_com = {}

    with open(structure_filename, "r") as structure_file:
        for line in structure_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                chain_id = line[20:23].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = atom_name[-1]

                mass = ATOMIC_MASSES.get(element, 0)
                if mass == 0:
                    continue

                if chain_id not in chain_com:
                    chain_com[chain_id] = {"mass": 0.0, "com": [0.0, 0.0, 0.0]}

                chain_com[chain_id]["mass"] += mass
                chain_com[chain_id]["com"][0] += x * mass
                chain_com[chain_id]["com"][1] += y * mass
                chain_com[chain_id]["com"][2] += z * mass

    for chain_id in chain_com:
        total_mass = chain_com[chain_id]["mass"]
        if total_mass > 0:
            chain_com[chain_id]["com"] = [c / total_mass for c in chain_com[chain_id]["com"]]

    return chain_com


def _normalize_element_name(element_name):
    if not element_name:
        return ""
    return str(element_name).strip().upper()


def _mass_from_element_name(element_name):
    normalized = _normalize_element_name(element_name)
    return ATOMIC_MASSES.get(normalized, 0)


def _altloc_allowed(atom):
    altloc = getattr(atom, "altloc", "")
    return not altloc or str(altloc) in {" ", "", "A"}


def calculate_center_of_mass_gemmi(structure_filename):
    if gemmi is None:
        raise RuntimeError("Gemmi is not installed. Install dependencies or use --parser legacy.")

    structure = gemmi.read_structure(structure_filename)
    chain_com = {}

    for model in structure:
        for chain in model:
            chain_id = str(chain.name).strip()
            if not chain_id:
                continue

            for residue in chain:
                for atom in residue:
                    if not _altloc_allowed(atom):
                        continue

                    element_name = _normalize_element_name(getattr(atom.element, "name", ""))
                    mass = _mass_from_element_name(element_name)
                    if mass == 0:
                        continue

                    if chain_id not in chain_com:
                        chain_com[chain_id] = {"mass": 0.0, "com": [0.0, 0.0, 0.0]}

                    chain_com[chain_id]["mass"] += mass
                    chain_com[chain_id]["com"][0] += atom.pos.x * mass
                    chain_com[chain_id]["com"][1] += atom.pos.y * mass
                    chain_com[chain_id]["com"][2] += atom.pos.z * mass

    for chain_id in chain_com:
        total_mass = chain_com[chain_id]["mass"]
        if total_mass > 0:
            chain_com[chain_id]["com"] = [c / total_mass for c in chain_com[chain_id]["com"]]

    return chain_com


def calculate_center_of_mass(structure_filename, parser="auto"):
    structure_path = Path(structure_filename)
    parser = parser.lower()

    if parser not in {"auto", "gemmi", "legacy"}:
        raise ValueError(f"Unsupported parser '{parser}'. Expected one of: auto, gemmi, legacy.")

    if parser == "legacy":
        print("Warning: using the legacy COM parser. It is a compatibility fallback and may be less reliable than Gemmi.")
        return calculate_center_of_mass_legacy(structure_filename)

    if parser == "gemmi":
        return calculate_center_of_mass_gemmi(structure_filename)

    if structure_path.suffix.lower() == ".cif":
        return calculate_center_of_mass_gemmi(structure_filename)

    if gemmi is not None:
        return calculate_center_of_mass_gemmi(structure_filename)

    print("Warning: Gemmi is unavailable; falling back to the legacy COM parser, which may be less reliable.")
    return calculate_center_of_mass_legacy(structure_filename)


def write_com_to_csv(chain_com, output_filename):
    output_path = Path(output_filename)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_filename, "w", newline="") as csvfile:
        fieldnames = ["chain", "x", "y", "z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for chain_id in sorted(chain_com):
            com_data = chain_com[chain_id]
            writer.writerow(
                {
                    "chain": chain_id,
                    "x": com_data["com"][0],
                    "y": com_data["com"][1],
                    "z": com_data["com"][2],
                }
            )

    print(f"Center of mass for each chain has been written to '{output_filename}'.")


def run(args):
    chain_com = calculate_center_of_mass(args.structure, parser=args.parser)
    write_com_to_csv(chain_com, args.output)


def build_parser(parser):
    parser.add_argument("--structure", required=True, help="Path to the input structure file.")
    parser.add_argument(
        "--parser",
        choices=["auto", "gemmi", "legacy"],
        default="auto",
        help="Parsing backend. 'auto' prefers Gemmi and falls back to the legacy compatibility parser for non-CIF input.",
    )
    parser.add_argument(
        "--output",
        default="center_of_mass.csv",
        help="Output CSV path for chain center-of-mass coordinates.",
    )
    parser.set_defaults(func=run)
    return parser


def main(argv=None):
    parser = argparse.ArgumentParser(description="Calculate center of mass for each chain.")
    build_parser(parser)
    args = parser.parse_args(argv)
    args.func(args)


if __name__ == "__main__":
    main()
