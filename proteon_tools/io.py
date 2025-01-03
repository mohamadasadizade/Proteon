"""
Module for Input/Output utilities related to protein structure analysis.
Includes functions for reading and writing PDB files, FASTA files, and other text-based formats.
"""

def readfile(filepath):
    """Reads a file and returns a list of lines.

    Args:
        filepath (str): The path to the file to be read.

    Returns:
        list: A list of lines from the file, each stripped of leading/trailing whitespaces.
    """
    with open(filepath, 'r') as f:
        return [line.strip() for line in f.readlines()]

def writefile(filepath, lines):
    """Writes a list of lines to a file.

    Args:
        filepath (str): The path to the file where the content will be written.
        lines (list): A list of lines to be written to the file.
    """
    with open(filepath, 'w') as f:
        for line in lines:
            f.write(line + '\n')

def pdb_to_list(filepath):
    """Reads a PDB file and returns a list of its lines.

    Args:
        filepath (str): The path to the PDB file.

    Returns:
        list: A list of PDB file lines, each stripped of leading/trailing whitespaces.
    """
    return readfile(filepath)

def fasta_dict(fasta_file):
    """Converts a FASTA formatted file into a dictionary.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are the FASTA labels and values are sequences.
    """
    with open(fasta_file, 'r') as f:
        FASTAdict = {}
        FASTAlabel = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                FASTAlabel = line
                FASTAdict[FASTAlabel] = ''
            else:
                FASTAdict[FASTAlabel] += line
    return FASTAdict


def split_pdb_by_chain(pdb_file, output_prefix):
    """Splits a multi-chain PDB file into separate PDB files for each chain.

    Args:
        pdb_file (str): The path to the input PDB file.
        output_prefix (str): The prefix to be used for the output files (e.g., 'output_chain_').
    """
    chains = {}

    # Read the PDB file and classify atoms by chain
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]  # The chain ID is found at the 22nd column
                if chain_id not in chains:
                    chains[chain_id] = []
                chains[chain_id].append(line)

    # Write separate files for each chain
    for chain_id, lines in chains.items():
        output_file = f"{output_prefix}_{chain_id}.pdb"
        with open(output_file, 'w') as file:
            for line in lines:
                file.write(line)
        print(f"Chain {chain_id} written to {output_file}")


def pdb_merge(file1, file2, output_file):
    """Merges two PDB files into a single output file.

    Args:
        file1 (str): The path to the first PDB file.
        file2 (str): The path to the second PDB file.
        output_file (str): The path to the output merged PDB file.
    """
    data1 = readfile(file1)
    data2 = readfile(file2)
    merged_data = data1 + [''] + data2
    writefile(output_file, merged_data)

def pdb_to_fasta(pdb_file, fasta_file):
    """Converts a PDB file to a FASTA sequence file.

    Args:
        pdb_file (str): The path to the input PDB file.
        fasta_file (str): The path to the output FASTA file.
    """
    sequence = ''

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                # Extract the chain ID and residue name (ignores heteroatoms like water)
                chain_id = line[21]
                residue_name = line[17:20].strip()
                # Include only standard amino acids
                if residue_name != 'HOH' and len(residue_name) == 3:
                    sequence += residue_name

    # Write the sequence to the FASTA file
    with open(fasta_file, 'w') as output_handle:
        output_handle.write(f">protein_sequence\n{sequence}")



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
