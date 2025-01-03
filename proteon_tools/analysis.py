"""
Include some functions for sequence-based protein analysis like residue composition. 
"""
import numpy as np
import mdtraj as md


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

fit_box = get_fit_box('native_34850_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb')
print(fit_box)




