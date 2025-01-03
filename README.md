# Proteon

`Proteon` is a Python package designed for protein sequence and structure analysis and manipulation. It includes utilities for protein sequence-based analysis, structure manipulation, and interaction with biological data formats like FASTA and PDB. This package aims to simplify the handling of protein sequences and structures, including residue composition analysis, protein property calculations, and PDB file manipulations.

## Features

### Sequence-Based Protein Analysis
- **Residue Composition**: Calculate residue composition from FASTA or PDB files.
- **Protein Properties**: Compute key protein properties such as molecular weight, net charge, and length.
- **Amino Acid Percentages**: Calculate the percentage of hydrophobic, polar, and charged residues in a protein sequence.

### Protein Structure Analysis
- **B-factor Plotting**: Visualize B-factors from PDB files to assess protein flexibility.
- **Fit Box Calculation**: Calculate an optimal fit box for a protein trajectory.

### Structure Manipulation
- **Chain Renumbering**: Change residue numbering in PDB files for specific chains.
- **Chain Renaming**: Rename chains in a PDB file.
- **PDB Merging**: Merge two PDB files into a single file.
- **Hydrogen Removal**: Remove hydrogen atoms from PDB files.

### Input/Output Utilities
- **FASTA to Dictionary**: Convert FASTA files to dictionaries.
- **PDB to FASTA**: Convert a PDB file to a FASTA sequence.
- **File Reading and Writing**: Simple utilities to read from and write to files.

## Installation

You can install `Proteon` from GitHub using `pip`:

```bash
git clone https://github.com/mohamadasadizade/Proteon.git
cd Proteon
pip install . 
