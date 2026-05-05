import csv

# Atomic masses (in atomic mass units, g/mol) for commonly used elements
atomic_masses = {
    'H': 1.008,
    'C': 12.011,
    'N': 14.007,
    'O': 15.999,
    'P': 30.974,
    'S': 32.06
    # Add more elements if needed
}

# Function to parse the PDB file and calculate the center of mass for each chain
def calculate_center_of_mass(pdb_filename):
    chain_com = {}

    # Open the PDB file for reading
    with open(pdb_filename, 'r') as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):  # Process only atoms
                # Extract atom details
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain_id = line[20:23].strip()  # Chain ID (this is the part we care about)
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                element = atom_name[-1]  # Assuming last character in atom name is element symbol
                
                # Look up the atomic mass
                mass = atomic_masses.get(element, 0)  # Default to 0 if element is not found

                # Skip atoms with no valid mass
                if mass == 0:
                    continue

                # Initialize chain data if not yet present
                if chain_id not in chain_com:
                    chain_com[chain_id] = {"mass": 0.0, "com": [0.0, 0.0, 0.0]}

                # Update the center of mass calculation for the chain
                chain_com[chain_id]["mass"] += mass
                chain_com[chain_id]["com"][0] += x * mass
                chain_com[chain_id]["com"][1] += y * mass
                chain_com[chain_id]["com"][2] += z * mass

    # Finalize the center of mass for each chain
    for chain_id in chain_com:
        total_mass = chain_com[chain_id]["mass"]
        if total_mass > 0:
            chain_com[chain_id]["com"] = [c / total_mass for c in chain_com[chain_id]["com"]]

    return chain_com

# Function to write the center of mass data to a CSV file
def write_com_to_csv(chain_com, output_filename="center_of_mass.csv"):
    with open(output_filename, "w", newline="") as csvfile:
        fieldnames = ["chain", "x", "y", "z"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for chain_id, com_data in chain_com.items():
            writer.writerow({
                "chain": chain_id,
                "x": com_data["com"][0],
                "y": com_data["com"][1],
                "z": com_data["com"][2]
            })

    print(f"Center of mass for each chain has been written to '{output_filename}'.")

# Example usage
pdb_filename = "asu_plus_surroundings_saved2icosbin2.pdb"  # Replace with your PDB file path
chain_com = calculate_center_of_mass(pdb_filename)
write_com_to_csv(chain_com)

