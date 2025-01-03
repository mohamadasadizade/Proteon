# Proteon

`Proteon` is a Python package for sequence-based analysis, structure manipulation, and interaction with biological data formats like FASTA and PDB. 

## Features

### Sequence-Based Protein Analysis
- **Residue Composition**: Calculate residue composition from FASTA or PDB files.
- **Protein Properties**: Compute key protein properties such as molecular weight, net charge, and length.
- **Amino Acid Percentages**: Calculate the percentage of hydrophobic, polar, and charged residues in a protein sequence.

### Protein Structure Analysis
- **B-factor Plotting**: Visualize B-factors from PDB files.
- **Fit Box Calculation**: An optimal fit box for simulation studies.

### PDB Manipulation
- **Chain Renumbering**: Change residue numbering in PDB files for specific chains.
- **Chain Renaming**: Rename chains in a PDB file.
- **PDB Merging**: Merge two PDB files into one file.
- **Hydrogen Removal**: Remove hydrogen atoms from PDB files.

### File handling
- **FASTA to Dictionary**: Convert FASTA files to dictionaries.
- **PDB to FASTA**: Convert a PDB file to a FASTA sequence.
- **File Reading and Writing**: Simple utilities to read from and write to files.

## Installation

You can install `Proteon` from GitHub using `pip`:

```bash
git clone https://github.com/mohamadasadizade/Proteon.git
cd Proteon
pip install . 
