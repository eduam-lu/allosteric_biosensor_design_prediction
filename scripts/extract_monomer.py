"""
Given a folder with multimeric protein structures in PDB or cif format, this script extracts each monomer to a another folder
"""
import pymol
from pymol import cmd
from pathlib import Path 
import argparse

parser = argparse.ArgumentParser(description="Extract monomers from multimeric PDB files using PyMOL.")
parser.add_argument("--input_folder", type=str, help="Path to the folder containing multimeric PDB files.")
parser.add_argument("--output_folder", type=str, help="Path to the folder to save extracted monomer PDB files.")
args = parser.parse_args()

Path(args.output_folder).mkdir(parents=True, exist_ok=True)
pymol.finish_launching(['pymol', '-qc'])  # Launch PyMOL in quiet mode

def extract_chains(pdb_path, output_path, chains_to_extract, label):
    """
    Use PyMOL to extract specified chains from a PDB and save to a new file.

    Parameters:
    - pdb_path: Path object
    - output_path: Path to folder
    - chains_to_extract: list of str (e.g., ['A', 'B'])
    - label: str, a label to append to the output filename (e.g., 'dimer')
    """
    base_name = pdb_path.stem
    output_file = Path(output_path) / f"{base_name}_{label}.pdb"

    cmd.load(str(pdb_path), 'structure')
    selection_str = " or ".join([f"chain {c}" for c in chains_to_extract])
    cmd.select("extracted", selection_str)
    cmd.save(str(output_file), "extracted")
    cmd.delete("all")

for CIF_file in Path(args.input_folder).glob("*.cif"):
    extract_chains(CIF_file, args.output_folder, ['A'], 'monomer')

for pdb_file in Path(args.input_folder).glob("*.pdb"):
    extract_chains(pdb_file, args.output_folder, ['A'], 'monomer')