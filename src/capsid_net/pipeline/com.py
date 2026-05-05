import argparse
import csv


ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.06,
}


def calculate_center_of_mass(structure_filename):
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


def write_com_to_csv(chain_com, output_filename):
    with open(output_filename, "w", newline="") as csvfile:
        fieldnames = ["chain", "x", "y", "z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for chain_id, com_data in chain_com.items():
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
    chain_com = calculate_center_of_mass(args.structure)
    write_com_to_csv(chain_com, args.output)


def build_parser(parser):
    parser.add_argument("--structure", required=True, help="Path to the input structure file.")
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
