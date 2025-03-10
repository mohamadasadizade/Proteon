"""
Include some functions for sequence-based protein analysis like residue composition. 
"""
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt


def fasta_rescompo(fasta_sequence):
    """ Calculate residue composition from a FASTA sequence."""
    residue_counts = {}
    for residue in fasta_sequence:
        residue_counts[residue] = residue_counts.get(residue, 0) + 1
    print(f"Residue Composition: {residue_counts}")
    return residue_counts


def pdb_rescompo(pdb_file):
    """
    Calculate the residue composition from a PDB file.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        dict: A dictionary with residues as keys and their counts as values.
    """
    residue_counts = {}

    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":  # C-alpha atoms indicate residues
                    residue = line[17:20].strip()  # Residue name
                    residue_counts[residue] = residue_counts.get(residue, 0) + 1
    except FileNotFoundError:
        print(f"Error: File '{pdb_file}' not found.")
        return {}
    except Exception as e:
        print(f"Error: {e}")
        return {}

    return residue_counts


def propro(sequence):
    """Summarize protein properties: weight, net charge, and length."""
    weights = {
        'A': 71.08, 'C': 103.14, 'D': 115.09, 'E': 129.12,
        'F': 147.18, 'G': 57.05, 'H': 137.14, 'I': 113.16,
        'K': 128.17, 'L': 113.16, 'M': 131.19, 'N': 114.11,
        'P': 97.12, 'Q': 128.13, 'R': 156.19, 'S': 87.08,
        'T': 101.11, 'V': 99.14, 'W': 186.21, 'Y': 163.18
    }
    charges = {
        'A': 0, 'C': 0, 'D': -1, 'E': -1, 'F': 0, 'G': 0, 'H': 0,
        'I': 0, 'K': +1, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0,
        'R': +1, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
    }

    weight = sum(weights.get(aa, 0) for aa in sequence)
    net_charge = sum(charges.get(aa, 0) for aa in sequence)
    length = len(sequence)

    print(f"Protein Properties - Weight: {weight:.2f}, Net Charge: {net_charge}, Length: {length}")
    return weight, net_charge, length


def aa_percentages(sequence):
    """Calculate the percentage of hydrophobic, polar, and charged residues."""
    hydrophobic = set("AILMFWYV")
    polar = set("STNQ")
    charged = set("DEKRH")

    total = len(sequence)
    if total == 0:
        print("Protein sequence is empty.")
        return 0, 0, 0

    hydrophobic_count = sum(1 for aa in sequence if aa in hydrophobic)
    polar_count = sum(1 for aa in sequence if aa in polar)
    charged_count = sum(1 for aa in sequence if aa in charged)

    hydrophobic_percentage = (hydrophobic_count / total) * 100
    polar_percentage = (polar_count / total) * 100
    charged_percentage = (charged_count / total) * 100

    print(f"Amino Acid Percentages - Hydrophobic: {hydrophobic_percentage:.2f}%, Polar: {polar_percentage:.2f}%, Charged: {charged_percentage:.2f}%")
    return hydrophobic_percentage, polar_percentage, charged_percentage


def get_fit_box(traj):
    # Load the trajectory using mdtraj
    traj = md.load(traj)
    # Get the minimum and maximum coordinates along the x, y, and z axes
    x_min, y_min, z_min = np.min(traj.xyz[0], axis=0)
    x_max, y_max, z_max = np.max(traj.xyz[0], axis=0)
    # Calculate the optimal box size
    box_x = x_max - x_min
    box_y = y_max - y_min
    box_z = z_max - z_min
    # Increase the box size by 1 nm from each side
    box_x += 2.0
    box_y += 2.0
    box_z += 2.0
    # Return the box vectors
    return np.array([[box_x, 0, 0], [0, box_y, 0], [0, 0, box_z]])



def parse_pdb(pdb_file):
    """Parse the PDB file to extract residue numbers and B-factors."""
    residue_numbers = []
    b_factors = []

    with open(pdb_file, 'r') as file:
        for line in file:
            # Look for lines starting with ATOM or HETATM (which contain atomic info)
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Residue number is at columns 23 to 26 (4 characters wide)
                # B-factor is at columns 61 to 66 (6 characters wide)
                residue_number = int(line[22:26].strip())
                b_factor = float(line[60:66].strip())
                
                # Append residue number and B-factor to lists
                residue_numbers.append(residue_number)
                b_factors.append(b_factor)

    return residue_numbers, b_factors

def plot_bfactor(pdb_file, output_file, cutoff, dpi=300):
    """Plot the B-factor from the PDB file with a cutoff line."""
    # Parse the PDB file to get residue numbers and B-factors
    residue_numbers, b_factors = parse_pdb(pdb_file)

    # Create a plot
    plt.plot(residue_numbers, b_factors, linestyle='-', color='blue')
    plt.axhline(y=cutoff, color='red', linestyle='--')  # Cutoff line

    # Highlight areas above the cutoff
    plt.fill_between(residue_numbers, b_factors, cutoff, where=(np.array(b_factors) >= cutoff), interpolate=True, color='lightsalmon')

    # Add labels and titles
    plt.xlabel('Residue Number')
    plt.ylabel('B-factor')
    plt.title(f'B-factor vs Residue Number for {pdb_file}')

    # Save the plot to a file
    plt.savefig(output_file, dpi=dpi)
    plt.show()



