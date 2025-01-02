"""
To manipulate PDB files, such as renumbering residues, changing chain names, and removing hydrogen atoms.
"""

def pdb_to_list(filepath):
    """Convert PDB lines to list elements."""
    with open(filepath) as f:
        return [i.strip() for i in f.readlines()]

def chain_renumbering(inputpath: str, outputPath: str, addnum: int, ch: str):
    """Change the numbering of atoms in a particular chain.

    Parameters:
        inputpath (str): The file path which should be renumbered.
        outputPath (str): The file path to save the final PDB file.
        addnum (int): The number to add to the residue numbers.
        ch (str): The name of the chain to renumber.
    """
    pdb_list = pdb_to_list(inputpath)
    with open(outputPath, 'w') as output:
        for line in pdb_list:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                if line[21] == ch:
                    new_line = line[:23] + str(int(line[23:26]) + addnum).rjust(3) + line[26:]
                    output.write(new_line + '\n')
                else:
                    output.write(line + '\n')
            else:
                output.write(line + '\n')

def pdb_chain_rename(inputpath: str, outputPath: str, start_line: int, last_line: int, ch_name: str):
    """Change chain name in defined lines.

    Parameters:
        inputpath (str): The directory of input file.
        outputPath (str): The directory of output file.
        start_line (int): The first line to change chain name.
        last_line (int): The last line to change chain name.
        ch_name (str): One letter name of the chosen chain.
    """
    pdb_list = pdb_to_list(inputpath)
    with open(outputPath, 'w') as output:
        for i, line in enumerate(pdb_list):
            if start_line - 1 <= i <= last_line - 1:
                new_line = line[:21] + ch_name + line[22:]
                output.write(new_line + '\n')
            else:
                output.write(line + '\n')

def pdb_merge(file1: str, file2: str, output_file: str):
    """Merge two PDB files into a single file."""
    with open(file1) as f1, open(file2) as f2:
        data1 = f1.read()
        data2 = f2.read()

    with open(output_file, 'w') as output:
        output.write(data1 + '\n' + data2)

def split_pdb_by_chain(pdb_file: str, output_prefix: str):
    """Split a multi-chain PDB file into separate PDB files for each chain.

    Parameters:
        pdb_file (str): The path to the input PDB file.
        output_prefix (str): The prefix to be used for the output files.
    """
    chains = {}

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]  # The chain ID is at position 22 (0-indexed)
                if chain_id not in chains:
                    chains[chain_id] = []
                chains[chain_id].append(line)

    for chain_id, lines in chains.items():
        output_file = f"{output_prefix}_{chain_id}.pdb"
        with open(output_file, 'w') as file:
            for line in lines:
                file.write(line)
        print(f"Chain {chain_id} written to {output_file}")

def pdb_remove_hydrogens(inputpath: str, outputPath: str):
    """Remove hydrogen atoms (ATOM/HETATM records where the atom name starts with 'H').

    Parameters:
        inputpath (str): The input PDB file path.
        outputPath (str): The output PDB file path to save without hydrogens.
    """
    pdb_list = pdb_to_list(inputpath)
    with open(outputPath, 'w') as output:
        for line in pdb_list:
            if not (line.startswith("ATOM") or line.startswith("HETATM")) or line[13:14] == "H":
                continue
            output.write(line + '\n')

def check_missing_residues(pdb_file: str):
    """Check for missing residues in the PDB file.

    Parameters:
        pdb_file (str): The input PDB file.
    """
    pdb_list = pdb_to_list(pdb_file)
    residues = {}
    missing_residues = []

    for line in pdb_list:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            chain_id = line[21]
            residue_id = int(line[22:26].strip())
            if (chain_id, residue_id) not in residues:
                residues[(chain_id, residue_id)] = line
            else:
                continue

    # Find missing residues (based on the chain and consecutive residue numbers)
    for (chain_id, residue_id) in sorted(residues.keys()):
        if (chain_id, residue_id + 1) not in residues:
            missing_residues.append((chain_id, residue_id + 1))

    return missing_residues
