"""
ligand_sampling.py --ligand --structure --output --help

DESCRIPTION

This script, given an input ligand and a structure with a different ligand bound to it, will sample different conformations and positions for that new ligand to replace the other.
This is a previous step to binding site redesign with RF diffussion AA/ Ligand MPNN

FLAGS

--ligand: SMILES format (prolly a file is better)

--structure: path to pdb or cif with ligand bound to it

--output: where the output folders should be generated

--help: extended help

WORKFLOW

A. Sample conformers with RD kit and save them as pdb files

B. Identify ligand positions

1. Manually through moving the ligand in pymol and saving the coordinates as a pdb
2. Not manually, through alignment with the previous ligand and stochastic tilting, rotating and inverting 

C. Prepare final pdb files

1. Load a position
2. Load conformer
3. Align conformer with the position (maybe find a way to specify atoms specially important to align)
4. Iterate

D. Final visual adjustments if needed

Eduardo Amo González
2025-2026
"""
### IMPORT MODULES ###########################################################################################################################

import argparse
import pandas as pd
from pathlib import Path
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import numpy as np
import pymol

### FUNCTIONS #################################################################################################################################

def sample_conformers_v0(molecule,n_conformers,output,ligand_name = "lig"):
    # Check if conformers where already sampled
    if Path(f"{output}/conformers/").exists():
        print("Conformers already found. Skipping conformer calculation")
        return
    # Load molecule as SMILE format
    mol = Chem.MolFromSmiles(molecule) # "c1ccccc1C(=O)O"
    mol = Chem.AddHs(mol)

    # Generate conformers
    params = AllChem.ETKDG()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)

    # Minimize each conformer
    for cid in conf_ids:
        AllChem.MMFFOptimizeMolecule(mol, confId=cid)

    # Save each conformer as a separate PDB
    Path(f"{output}/conformers/").mkdir(parents=True, exist_ok=True)

    for cid in conf_ids:
        filename = f"{output}/conformers/{ligand_name}_conf_{cid}.pdb"
        Chem.MolToPDBFile(mol, filename, confId=cid)
    print (f"{n_conformers} were saved succesfully")
    return

def sample_conformers(molecule, n_conformers, rmsd_cutoff, output, ligand_name="lig"):
    out_dir = Path(f"{output}/conformers/")
    
    # Check if conformers already exist
    if out_dir.exists():
        print("Conformers already found. Skipping conformer calculation.")
        return

    # Prepare output folder
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load molecule
    mol = Chem.MolFromSmiles(molecule)
    mol = Chem.AddHs(mol)

    # Generate conformers
    print("Generating conformers...")
    params = AllChem.ETKDG()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)

    energies = []
    print("Minimizing conformers...")

    # Minimize each conformer & store energies
    for cid in conf_ids:
        prop = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        ff = AllChem.MMFFGetMoleculeForceField(mol, prop, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        energies.append(energy)

    energies = np.array(energies)

    # Save PDBs
    for cid in conf_ids:
        filename = out_dir / f"{ligand_name}_conf_{cid}.pdb"
        Chem.MolToPDBFile(mol, str(filename), confId=cid)

    print(f"{n_conformers} conformers generated and saved successfully.")

    # STATISTICS

    # Pairwise RMSD matrix
    rmsd_matrix = np.zeros((n_conformers, n_conformers))
    for a, i in enumerate(conf_ids):
        for b, j in enumerate(conf_ids):
            if b <= a: 
                continue
            rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=i, refId=j)
            rmsd_matrix[a, b] = rmsd
            rmsd_matrix[b, a] = rmsd

    # RMSD to lowest-energy conformer
    best_idx = energies.argmin()
    best_conf_id = conf_ids[int(best_idx)]
    lowest_energy_path = out_dir / f"{ligand_name}_conf_{best_conf_id}.pdb"

    rmsd_to_best = np.array([
        rdMolAlign.GetBestRMS(mol, mol, prbId=conf_ids[i], refId=best_conf_id)
        for i in range(n_conformers)
    ])

    # Count unique conformers using RMSD threshold
    
    unique = [0]  # always include lowest-energy conformer
    for i in range(1, n_conformers):
        if all(rmsd_matrix[i, j] > rmsd_cutoff for j in unique):
            unique.append(i)

    # Statistics output
    stats_df = pd.DataFrame({
        "conf_id": conf_ids,
        "energy": energies,
        "rmsd_to_lowest": rmsd_to_best,
        "is_unique": [cid in unique for cid in conf_ids]
    })

    stats_df.to_csv(out_dir / "conformer_stats.csv", index=False)
    np.savetxt(out_dir / "rmsd_matrix.txt", rmsd_matrix, fmt="%.3f")
    print(f"Unique conformers (RMSD > {rmsd_cutoff} Å): {len(unique)}")

    return lowest_energy_path

def identify_ligand_positions(molecule_path, structure, num_positions, output):

    # Check if the user has done it manually
    out_dir = Path(f"{output}/positions/")
    
    # Check if conformers already exist
    if out_dir.exists():
        print("Positions already found. Skipping position calculation.")
        return
    
    # Identify ligand positions through alignment and stochastic transformations
    # A. Load the new ligand in the structure
    
    # B. Align the new ligand with the old one

    # C. Remove old ligand

    # D. Apply stochastic tilting, rotating and inverting

    # E. Save positions

    return



### PARAMS ####################################################################################################################################

# Conformer sampling parameters
num_conformers = 5
conformer_rmsd_cutoff = 0.75

# Position sampling parameters
num_positions = 5
### INPUT CHECK ###############################################################################################################################

parser = argparse.ArgumentParser(
    description="This script runs the pipeline shapedesign with an additional RF_diffusion step"
)
parser.add_argument('--ligand', help="SMILES code for the ligand", type=str)
parser.add_argument('--structure', help="Path to the structure to redesign", type=str)
parser.add_argument('--output', help="Folder where the outputs will be stored", type=str)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# If --ligand, structure or output was not provided, show error and exit
if not args.ligand:
    parser.error("--ligand is required unless --detailed-help is used")

if not args.structure:
    parser.error("--structure is required unless --detailed-help is used")

if not args.output:
    parser.error("--output is required unless --detailed-help is used")

# Assign variables
ligand_molecule = str(args.ligand)
structure_path = Path(args.structure)
output_path = Path(args.output)

# If structure not found, error and exit
if not structure_path.exists():
    parser.error(f"The file '{structure_path}' does not exist.")

if structure_path.is_dir():
    parser.error(f"The path '{structure_path}' is a directory. File required")

# Create output folder if it doesn't exist

if not output_path.exists():
    output_path.mkdir(parents=True, exist_ok=True)
else:
    if not output_path.is_dir():
        parser.error(f"The path '{output_path}' is not a directory")

### MAIN EXECUTION ############################################################################################################################

### A. Sample conformers with RDkit
lowest_energy_conformer = sample_conformers(ligand_molecule,num_conformers,conformer_rmsd_cutoff,output_path,"Tc")


### B. Identify ligand positions
identify_ligand_positions(lowest_energy_conformer,structure_path,num_positions,output_path)

### C. Prepare final pdb files
