"""
Given a folder full of PDBs or cifs from a cofolding prediction, it checks how many of them succeeded

To do this, it is given a list of residues and the ligand name. It computes the centre of mass of both the group and the ligand, computes the distance 
between them, and checks if it is below a certain threshold (e.g. 5 Å). If it is, it counts as a success, otherwise as a failure. 
It then prints the number of successes and failures, and the success rate.

Also, results should be saved in a csv file with columns: file name, success (True/False), distance (float).

python check_cofolding_success.py \
    --input_dir ./predictions \
    --residues A:10 A:25 \
    --ligand LIG \
    --threshold 4.5 \
    --output success_metrics.csv
"""
import os
import glob
import argparse
import numpy as np
import pandas as pd
import warnings
from Bio.PDB import PDBParser, MMCIFParser
from Bio import BiopythonWarning

# Suppress minor parsing warnings from Biopython
warnings.simplefilter('ignore', BiopythonWarning)

def parse_residue_list(res_list_str):
    """
    Parses a list of residues in the format Chain:ResNum (e.g., A:10 B:25)
    Returns a list of tuples: [('A', 10), ('B', 25)]
    """
    target_residues = []
    for res in res_list_str:
        try:
            chain, num = res.split(':')
            target_residues.append((chain, int(num)))
        except ValueError:
            print(f"Warning: Ignored invalid residue format '{res}'. Use 'Chain:Number' (e.g., A:15).")
    return target_residues

def get_com(atoms):
    """Calculate the Center of Mass for a given list of Bio.PDB atoms."""
    if not atoms:
        return None
    coords = [atom.coord for atom in atoms]
    return np.mean(coords, axis=0)

def evaluate_structure(file_path, target_residues, ligand_name):
    """
    Parses the structure, calculates the distance between the COM of the 
    specified residues and the COM of the ligand.
    """
    # Initialize appropriate parser
    if file_path.endswith('.pdb'):
        parser = PDBParser(QUIET=True)
    elif file_path.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        return None
        
    structure = parser.get_structure('struct', file_path)
    
    res_atoms = []
    ligand_atoms = []
    
    # We only look at the first model (model 0)
    for model in structure:
        for chain in model:
            for res in chain:
                # Check if it's one of the target residues
                if (chain.id, res.id[1]) in target_residues:
                    res_atoms.extend(list(res.get_atoms()))
                
                # Check if it's the target ligand
                # strip() is used because PDB format strings sometimes have padding
                if res.resname.strip() == ligand_name:
                    ligand_atoms.extend(list(res.get_atoms()))
        break 

    res_com = get_com(res_atoms)
    lig_com = get_com(ligand_atoms)
    
    if res_com is None or lig_com is None:
        return np.nan # Cannot compute distance if entities are missing
        
    # Calculate Euclidean distance between the two Centers of Mass
    distance = np.linalg.norm(res_com - lig_com)
    return distance

def main():
    parser = argparse.ArgumentParser(description="Evaluate cofolding success by measuring distance between residue COM and Ligand COM.")
    parser.add_argument("-i", "--input_dir", required=True, help="Path to the folder containing PDB/CIF files.")
    parser.add_argument("-r", "--residues", required=True, nargs='+', help="List of residues (Format: CHAIN:RESNUM, e.g., A:10 B:45).")
    parser.add_argument("-l", "--ligand", required=True, help="Name of the ligand (e.g., LIG, ATP).")
    parser.add_argument("-t", "--threshold", type=float, default=5.0, help="Distance threshold in Ångströms for success (default: 5.0).")
    parser.add_argument("-o", "--output", default="cofolding_results.csv", help="Output CSV file name (default: cofolding_results.csv).")
    
    args = parser.parse_args()
    
    target_residues = parse_residue_list(args.residues)
    if not target_residues:
        print("Error: No valid target residues provided. Exiting.")
        return

    # Find all PDB and CIF files in the target directory
    search_pdb = os.path.join(args.input_dir, "*.pdb")
    search_cif = os.path.join(args.input_dir, "*.cif")
    files = glob.glob(search_pdb) + glob.glob(search_cif)
    
    if not files:
        print(f"No .pdb or .cif files found in directory: {args.input_dir}")
        return

    print(f"Found {len(files)} structure files. Processing...")

    results = []
    successes = 0
    failures = 0

    for file_path in files:
        file_name = os.path.basename(file_path)
        distance = evaluate_structure(file_path, target_residues, args.ligand)
        
        # If distance is NaN (ligand or residues were missing), it automatically fails
        if pd.isna(distance):
            success = False
            failures += 1
        else:
            success = bool(distance <= args.threshold)
            if success:
                successes += 1
            else:
                failures += 1
                
        results.append({
            'file name': file_name,
            'success': success,
            'distance': round(float(distance), 3) if not pd.isna(distance) else distance
        })

    # Print summary
    total = successes + failures
    success_rate = (successes / total) * 100 if total > 0 else 0.0
    
    print("\n--- Summary ---")
    print(f"Total Processed : {total}")
    print(f"Successes       : {successes}")
    print(f"Failures        : {failures}")
    print(f"Success Rate    : {success_rate:.2f}%")
    print("---------------")

    # Save to CSV
    df = pd.DataFrame(results)
    df.to_csv(args.output, index=False)
    print(f"\nResults saved to '{args.output}'")

if __name__ == "__main__":
    main()