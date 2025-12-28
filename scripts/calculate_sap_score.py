"""
calculate_sap_score --input_folder <input_folder> --output_folder <output_folder>

Given a folder with PDB files, returns a csv with the SAP scores for each structure.
SAP (Spatial Aggregation Propensity) is calculated based on the hydrophobicity and solvent accessibility
of residues within a certain radius.
"""

### MODULES ####################################################################
import os
import glob
import pandas as pd
import numpy as np
import Bio.PDB
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.NeighborSearch import NeighborSearch
import argparse
from pathlib import Path
### PARAMETERS & CONSTANTS #####################################################
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

# Black & Mould Hydrophobicity Scale (Normalized)
# Values range from -1 (Hydrophilic) to 1 (Hydrophobic)
# Glycine is usually 0.
HYDROPHOBICITY = {
    'PHE': 1.000, 'LEU': 0.943, 'ILE': 0.943, 'TYR': 0.880,
    'TRP': 0.849, 'VAL': 0.825, 'MET': 0.738, 'CYS': 0.680,
    'ALA': 0.500, 'THR': 0.450, 'SER': 0.300, 'GLY': 0.000,
    'HIS': -0.100, 'LYS': -0.200, 'ARG': -0.300, 'GLN': -0.400,
    'ASN': -0.500, 'GLU': -0.600, 'ASP': -0.700, 'PRO': 0.000 # Proline often treated as neutral/structural
}

# Max Sidechain SASA (Approximate values from Ala-X-Ala tripeptides)
# Used to normalize the actual SASA to get a fraction (0.0 to 1.0)
MAX_SASA_SIDECHAIN = {
    'ALA': 67.0,  'ARG': 196.0, 'ASN': 113.0, 'ASP': 106.0,
    'CYS': 104.0, 'GLN': 144.0, 'GLU': 138.0, 'GLY': 0.1, # Gly has no sidechain, strictly speaking
    'HIS': 151.0, 'ILE': 140.0, 'LEU': 137.0, 'LYS': 167.0,
    'MET': 160.0, 'PHE': 175.0, 'PRO': 105.0, 'SER': 80.0,
    'THR': 102.0, 'TRP': 217.0, 'TYR': 187.0, 'VAL': 117.0
}

radius = 5.0  # Radius in Angstroms for neighbor search. Default is 5.0 Å since it's what Rosetta uses.

### FUNCTIONS ######################################################################

def get_sidechain_sasa(residue):
    """Sum SASA of all atoms in residue excluding backbone (N, C, O, CA)."""
    sasa = 0.0
    for atom in residue:
        if atom.name not in ['N', 'C', 'O', 'CA']:
            sasa += getattr(atom, 'sasa', 0.0)
    return sasa

def process_single_pdb_rosetta_style(pdb_path):
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_path)
    except:
        return None

    # 1. Calculate SASA for all atoms
    sr = ShrakeRupley()
    sr.compute(structure, level="A")

    # 2. Prepare Data Structures
    all_sidechain_atoms = []
    atom_potentials = {} 
    residue_map = {}     

    for model in structure:
        for chain in model:
            for residue in chain:
                if not Bio.PDB.is_aa(residue): continue
                
                res_name = residue.get_resname()
                hydro = HYDROPHOBICITY.get(res_name, 0.0)
                max_sasa = MAX_SASA_SIDECHAIN.get(res_name, 100.0)
                
                scale_factor = hydro / max_sasa if max_sasa > 0.1 else 0.0

                current_res_atoms = []
                
                for atom in residue:
                    if atom.name in ['N', 'C', 'O', 'OXT']: continue
                    
                    all_sidechain_atoms.append(atom)
                    current_res_atoms.append(atom)
                    
                    atom_sasa = getattr(atom, 'sasa', 0.0)
                    atom_potentials[atom] = atom_sasa * scale_factor

                residue_map[residue] = current_res_atoms

    # 3. Build KDTree (radius = 5.0 matches Rosetta)
    if not all_sidechain_atoms: return None
    ns = NeighborSearch(all_sidechain_atoms)

    sap_scores = []
    sequence_list = []

    # 4. Calculate SAP per Residue
    for model in structure:
        for chain in model:
            for residue in chain:
                if not Bio.PDB.is_aa(residue): continue

                res_name_3 = residue.get_resname()
                sequence_list.append(AA_3TO1.get(res_name_3, 'X'))

                res_sap_sum = 0.0
                my_atoms = residue_map.get(residue, [])
                
                for my_atom in my_atoms:
                    neighbors = ns.search(my_atom.get_coord(), 5.0, level='A')
                    for neighbor_atom in neighbors:
                        res_sap_sum += atom_potentials.get(neighbor_atom, 0.0)

                sap_scores.append(res_sap_sum)

    # 5. Format Output
    scores = np.array(sap_scores)
    file_id = os.path.basename(pdb_path).replace(".pdb", "").split("_seq_")[0]
    
    sap_string = ",".join([f"{s:.3f}" for s in scores])
    seq_string = "".join(sequence_list)

    return {
        "pdb_id": file_id,
        "max_sap": np.max(scores),
        "min_sap": np.min(scores),                 # <--- Added
        "sum_positive_sap": np.sum(scores[scores > 0]),
        "global_sap_sum": np.sum(scores),
        "avg_sap": np.mean(scores),
        "residue_count": len(scores),              # <--- Added
        "sequence": seq_string,
        "sap_string": sap_string
    }

def process_single_pdb(pdb_path, radius=10.0):
    """
    Calculates SAP metrics for a single PDB file.
    Returns a dictionary of metrics including the per-residue SAP string.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_path)
    except Exception as e:
        print(f"Error parsing {pdb_path}: {e}")
        return None

    # 1. Calculate Surface Area (SASA)
    sr = ShrakeRupley()
    sr.compute(structure, level="A")

    # 2. Collect residues and CA atoms
    atoms_for_search = []
    residues = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                if Bio.PDB.is_aa(residue) and 'CA' in residue:
                    residues.append(residue)
                    atoms_for_search.append(residue['CA'])

    if not residues:
        return None

    ns = NeighborSearch(atoms_for_search)
    sap_scores = []
    sequence_list = [] # To store the 1-letter sequence

    # 3. Calculate SAP for each residue
    for target_res in residues:
        # -- Sequence Handling --
        res_name_3 = target_res.get_resname()
        res_name_1 = AA_3TO1.get(res_name_3, 'X') # Default to X if unknown
        sequence_list.append(res_name_1)
        
        # -- SAP Calculation --
        target_sap = 0.0
        neighbors = ns.search(target_res['CA'].get_coord(), radius, level='R')
        
        for neighbor in neighbors:
            n_name = neighbor.get_resname()
            
            # Calculate Fraction Exposed
            actual_area = get_sidechain_sasa(neighbor)
            max_area = MAX_SASA_SIDECHAIN.get(n_name, 100.0)
            
            perc_exposed = 0.0
            if max_area >= 1.0:
                perc_exposed = min(actual_area / max_area, 1.0)
            
            hydrophobicity = HYDROPHOBICITY.get(n_name, 0.0)
            target_sap += (perc_exposed * hydrophobicity)
        
        sap_scores.append(target_sap)

    # 4. Compile Metrics
    scores = np.array(sap_scores)
    file_id = os.path.basename(pdb_path).replace(".pdb", "")
    if "_seq_" in file_id:
        file_id = file_id.split("_seq_")[0]

    # Create the SAP String (comma separated values, rounded to 3 decimals)
    # Example: "0.051,1.200,-0.342,..."
    sap_string_formatted = ",".join([f"{s:.3f}" for s in scores])
    
    # Create the Sequence String
    # Example: "MKTIIAL..."
    sequence_string = "".join(sequence_list)

    return {
        "pdb_id": file_id,
        "global_sap_sum": np.sum(scores),
        "sum_positive_sap": np.sum(scores[scores > 0]),
        "max_sap": np.max(scores),
        "min_sap": np.min(scores),
        "avg_sap": np.mean(scores),
        "residue_count": len(scores),
        "sequence": sequence_string,    # The AA sequence
        "sap_string": sap_string_formatted # The scores matching the sequence
    }


def analyze_folder(folder_path, output_csv="sap_results.csv", radius=radius):
    """
    Iterates through a folder of PDBs, calculates SAP metrics,
    and returns a Pandas DataFrame.
    """
    pdb_files = glob.glob(os.path.join(folder_path, "*.pdb"))
    print(f"Found {len(pdb_files)} PDB files in {folder_path}...")
    
    results = []
    
    for i, pdb_file in enumerate(pdb_files):
        # Optional: Print progress every 10 files
        if (i+1) % 10 == 0:
            print(f"Processing {i+1}/{len(pdb_files)}...")
            
        metrics = process_single_pdb_rosetta_style(pdb_file)
        if metrics:
            results.append(metrics)
    
    # Create DataFrame
    df = pd.DataFrame(results)
    
    # Reorder columns for readability (put ID first)
    cols = ["pdb_id", "max_sap", "sum_positive_sap", "global_sap_sum", "avg_sap", "min_sap", "sap_string", "sequence", "residue_count"]
    df = df[cols]
    
    # Save to CSV
    if output_csv:
        df.to_csv(output_csv, index=False)
        print(f"Results saved to {output_csv}")
        
    return df


### INPUT CHECK #######################################################################

parser = argparse.ArgumentParser(description="Calculate SAP scores for PDB files in a folder.")
parser.add_argument('--input_folder', type=str, required=True, help='Path to folder containing PDB files.')
parser.add_argument('--output_folder', type=str, default='.', help='Path to folder to save output CSV.')    
args = parser.parse_args()

Path(args.output_folder).mkdir(parents=True, exist_ok=True)

### EXECUTION #########################################################################
if __name__ == "__main__":
    
    # Run analysis
    df_results = analyze_folder(args.input_folder, output_csv=os.path.join(args.output_folder, "sap_results.csv"))
    
    # Display the first few rows
    print("\nDataFrame Head:")
    print(df_results.head())