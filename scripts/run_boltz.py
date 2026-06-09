"""
Given a sequence fasta and a SMILES string, this script runs Boltz on the sequence and ligand. You can choose to run several iterations as well.

usage: python run_boltz.py --fasta_path <path_to_fasta> --smiles_string <smiles_string> --output_dir <output_directory> 
    --run_gnina --num_iterations <number_of_iterations> --ligand_name <ligand_name_for_output_files> --pocket_residues <list_of_pocket_residues>
"""

### IMPORTS
import argparse
import re
import subprocess
import os
import yaml
import json
import numpy as np
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import Bio.PDB as PDB
from Bio.PDB import PDBParser, MMCIFParser, DSSP, NeighborSearch, Selection
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.SASA import ShrakeRupley
import biotite.structure.io as bsio

### PARAMS
# Boltz2 Parameters
BOLTZ_DEVICES = 1
BOLTZ_RECYCLING_STEPS = 3
BOLTZ_SAMPLING_STEPS = 100
BOLTZ_DIFFUSION_SAMPLES = 1
BOLTZ_OUTPUT_FORMAT = 'pdb'
BOLTZ_SAMPLING_STEPS_AFFINITY = 100
BOLTZ_USE_MSA = True
BOLTZ_USE_FORCES = False
BOLTZ_NO_KERNELS = True
BOLTZ_POCKET_MAX_DIST = 5.0

# Gnina Parameters
GNINA_PATH = "/proj/berzelius-2023-361/users/x_eduam/gnina/build/bin/gnina"  # Path to gnina executable
GNINA_EXHAUSTIVENESS = 16
GNINA_BOX_SIZE = 15
GNINA_CNN_MODEL = "crossdock_default2018"  # or other gnina CNN model names

# Ligand Sampling Parameters
N_CONFORMERS = 10  # Number of conformers to generate
RMSD_CUTOFF = 0.5  # RMSD cutoff for conformer clustering
LIGAND_NAME = "LIG"  # Name of the ligand for output files

# Pocket Parameters
POCKET_RESIDUES = []  # List of pocket residues: [[chain, res_idx], ...], empty list for no pocket constraint

# Path Parameters
BOLTZ_CONDA_ENV = "boltz"

### FUNCTIONS


def boltz_yaml_generator_w_msa(row, yaml_path, ligand_smiles, pocket_list, max_dist=5.0):
    """
    Generates a Boltz-compatible YAML file for protein-ligand complexes.
    
    Args:
        row (dict/Series): Must contain 'file_ID' and 'sequence'.
        yaml_path (str): Directory to save the YAML.
        ligand_smiles (str): SMILES string for the ligand.
        pocket_list (list): List of contacts [[CHAIN, RES_IDX], ...]
        max_dist (float): Maximum distance for pocket constraints.
    """
    file_id = row['file_ID']
    seq = row['sequence']
    
    # Constructing the dictionary structure
    if pocket_list:
        yaml_data = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": seq 
                    }
                },
                {
                    "ligand": {
                        "id": "B",
                        "smiles": ligand_smiles
                    }
                }
            ],
            "constraints": [
                {
                    "pocket": {
                        "binder": "A",  # Defining the pocket on the Protein (A)
                        "contacts": pocket_list,
                        "max_distance": max_dist,
                        "force": False
                    }
                }
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "B"  # Calculate affinity for the Ligand (B)
                    }
                }
            ]
        }
    else:
        yaml_data = {
            "version": 1,
            "sequences": [
                {
                    "protein": {
                        "id": "A",
                        "sequence": seq
                    }
                },
                {
                    "ligand": {
                        "id": "B",
                        "smiles": ligand_smiles
                    }
                }
            ],
            "properties": [
                {
                    "affinity": {
                        "binder": "B"  # Calculate affinity for the Ligand (B)
                    }
                }
            ]
        }

    # Ensure output directory exists
    os.makedirs(yaml_path, exist_ok=True)
    file_output = os.path.join(yaml_path, f"{file_id}.yaml")

    with open(file_output, "w") as f:
        # sort_keys=False preserves the order of the dictionary
        yaml.dump(yaml_data, f, default_flow_style=False, sort_keys=False)
    
    return file_output

def run_Boltz2(row, ligand_smiles, output_path, pocket_list=None, max_dist=None, 
                use_msa=None, use_forces=None, no_kernels=None,
                devices=None, recycling_steps=None, sampling_steps=None, 
                diffusion_samples=None, output_format=None, sampling_steps_affinity=None):
    # Use global params if not provided
    if pocket_list is None:
        pocket_list = []
    if max_dist is None:
        max_dist = BOLTZ_POCKET_MAX_DIST
    if use_msa is None:
        use_msa = BOLTZ_USE_MSA
    if use_forces is None:
        use_forces = BOLTZ_USE_FORCES
    if no_kernels is None:
        no_kernels = BOLTZ_NO_KERNELS
    if devices is None:
        devices = BOLTZ_DEVICES
    if recycling_steps is None:
        recycling_steps = BOLTZ_RECYCLING_STEPS
    if sampling_steps is None:
        sampling_steps = BOLTZ_SAMPLING_STEPS
    if diffusion_samples is None:
        diffusion_samples = BOLTZ_DIFFUSION_SAMPLES
    if output_format is None:
        output_format = BOLTZ_OUTPUT_FORMAT
    if sampling_steps_affinity is None:
        sampling_steps_affinity = BOLTZ_SAMPLING_STEPS_AFFINITY
    
    # Generate the yaml for this row
    yaml_path = boltz_yaml_generator_w_msa(row, output_path, ligand_smiles, pocket_list, max_dist)
    # Base command with mandatory flags
    command_parts = [
        f"conda run -n {BOLTZ_CONDA_ENV} boltz predict",
        f"{yaml_path}",
        f"--out_dir={output_path}",
        f"--devices={devices}",
        f"--recycling_steps={recycling_steps}",
        f"--sampling_steps={sampling_steps}",
        f"--diffusion_samples={diffusion_samples}",
        f"--output_format={output_format}",
        f"--sampling_steps_affinity={sampling_steps_affinity}"
    ]

    # Conditional flags
    if use_msa:
        command_parts.append("--use_msa_server")
    
    if use_forces:
        command_parts.append("--use_potentials")
        
    if no_kernels:
        command_parts.append("--no_kernels")

    # Join all parts into a single string
    boltz_command = " ".join(command_parts)
    
    # Run the process
    boltz_process = subprocess.run(boltz_command, shell=True)
    # Remove yaml file

    return boltz_process

def parse_fasta(fasta_path):
    """
    Parses a fasta file and returns a list of dictionaries with 'header' and 'sequence'.
    
    Args:
        fasta_path (str): Path to the input fasta file.
    Returns:
        List[Dict]: A list of dictionaries, each containing 'header' and 'sequence'.
    """
    sequences = []
    with open(fasta_path, "r") as f:
        file_id = None
        sequence = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if file_id is not None:
                    sequences.append({"header": file_id, "sequence": sequence})
                file_id = line[1:]  # Remove the '>' character
                sequence = ""  # Reset sequence for the new entry
            else:
                sequence += line  # Append to the current sequence
        # Add the last entry after the loop ends
        if file_id is not None:
            sequences.append({"header": file_id, "sequence": sequence})
    return sequences

def extract_perresidue_plddt(file_path):
    """
    Extracts per-residue (per-position) pLDDT scores from a PDB file.
    
    In PDB files, B-factors are stored per-atom. This function extracts the 
    CA (alpha carbon) B-factor for each residue, which represents the pLDDT
    confidence score for that residue position.
    
    Returns:
        list[float]: List of pLDDT scores, one per residue position.
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('struct', file_path)
    except Exception as e:
        print(f"Error parsing PDB {file_path}: {e}")
        return []
    
    plddt_values = []
    
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip non-standard residues (waters, ligands, etc.)
                if not PDB.is_aa(residue):
                    continue
                
                # Extract CA (alpha carbon) B-factor for this residue
                if 'CA' in residue:
                    ca_atom = residue['CA']
                    b_factor = ca_atom.get_bfactor()
                    plddt_values.append(float(b_factor))
                else:
                    # If CA is missing, try to use average of all atoms
                    atom_bfactors = []
                    for atom in residue:
                        atom_bfactors.append(atom.get_bfactor())
                    
                    if atom_bfactors:
                        avg_bfactor = np.mean(atom_bfactors)
                        plddt_values.append(float(avg_bfactor))
    
    return plddt_values

def process_boltz_output(output_path, pdb_search_path=None):
    """
    Aggregate Boltz output metrics from a Boltz prediction folder.

    Args:
        output_path (str): Root Boltz output folder where confidence and affinity JSON files are stored.
        pdb_search_path (str, optional): Root folder for generated PDB files. If omitted it uses output_path.

    Returns:
        dict: {
            'avg_confidence_score': float | None,
            'avg_affinity_pred_value': float | None,
            'avg_affinity_probability_binary': float | None,
            'avg_plddt_per_position': list[float|None],
            'num_confidence_files': int,
            'num_affinity_files': int,
            'num_pdbs': int
        }
    """
    root_path = Path(output_path)
    search_root = Path(pdb_search_path) if pdb_search_path else root_path

    confidence_scores = []
    affinity_values = []
    affinity_binaries = []
    affinity_files = 0

    for conf_file in root_path.rglob('confidence_*.json'):
        try:
            with open(conf_file, 'r') as f:
                data = json.load(f)
        except Exception:
            continue

        confidence = data.get('confidence_score', data.get('overall_confidence', None))
        if confidence is not None:
            try:
                confidence_scores.append(float(confidence))
            except (TypeError, ValueError):
                pass

        file_id = conf_file.stem.replace('confidence_', '').replace('_model_0', '')
        affinity_file = conf_file.parent / f'affinity_{file_id}.json'
        if affinity_file.exists():
            affinity_files += 1
            try:
                with open(affinity_file, 'r') as af:
                    affinity_data = json.load(af)
            except Exception:
                continue

            pred_val = affinity_data.get('affinity_pred_value', affinity_data.get('affinity_pred', None))
            binary_val = affinity_data.get('affinity_probability_binary', affinity_data.get('affinity_binary', affinity_data.get('affinity_probability', None)))

            if pred_val is not None:
                try:
                    affinity_values.append(float(pred_val))
                except (TypeError, ValueError):
                    pass

            if binary_val is not None:
                try:
                    affinity_binaries.append(float(binary_val))
                except (TypeError, ValueError):
                    pass

    plddt_lists = []
    pdb_files = sorted(search_root.rglob('*.pdb'))
    for pdb_path in pdb_files:
        plddt_values = extract_perresidue_plddt(str(pdb_path))
        if isinstance(plddt_values, (list, tuple)) and plddt_values:
            try:
                arr = np.array(plddt_values, dtype=float)
                if arr.size > 0:
                    plddt_lists.append(arr)
            except Exception:
                continue

    avg_plddt_per_position = []
    if plddt_lists:
        max_len = max(arr.shape[0] for arr in plddt_lists)
        stacked = np.full((len(plddt_lists), max_len), np.nan, dtype=float)
        for idx, arr in enumerate(plddt_lists):
            stacked[idx, : arr.shape[0]] = arr

        mean_values = np.nanmean(stacked, axis=0)
        avg_plddt_per_position = [float(v) if not np.isnan(v) else None for v in mean_values]

    avg_confidence_score = float(np.mean(confidence_scores)) if confidence_scores else None
    avg_affinity_pred_value = float(np.mean(affinity_values)) if affinity_values else None
    avg_affinity_probability_binary = float(np.mean(affinity_binaries)) if affinity_binaries else None

    return {
        'avg_confidence_score': avg_confidence_score,
        'avg_affinity_pred_value': avg_affinity_pred_value,
        'avg_affinity_probability_binary': avg_affinity_probability_binary,
        'avg_plddt_per_position': avg_plddt_per_position,
        'num_confidence_files': len(confidence_scores),
        'num_affinity_files': affinity_files,
        'num_pdbs': len(pdb_files)
    }

def gnina_box_generator(pdb_file, residues, size=None):
    """
    Generates a box dictionary for gnina based on the provided residues.
    
    :param pdb_file: Path to the PDB file.
    :param residues: List of tuples [(chain, res_idx), ...] defining the pocket residues.
    :param size: Size of the box in Angstroms (default from GNINA_BOX_SIZE).
    :return: Dictionary with box parameters for gnina.
    """
    if size is None:
        size = GNINA_BOX_SIZE
    
    # Load structure and filter for CA atoms
    struct = bsio.load_structure(pdb_file)
    ca = struct[struct.atom_name == "CA"]
    
    # Extract coordinates of specified residues
    coords = []
    for residue in residues:
        if isinstance(residue, int):
            chain, res_idx = "A", residue
        else:
            chain, res_idx = residue
        mask = (ca.chain_id == chain) & (ca.res_id == res_idx)
        if np.any(mask):
            coords.append(ca.coord[mask][0])  # Assuming one CA per residue
    
    if not coords:
        raise ValueError("No valid residues found for box generation.")
    
    coords = np.array(coords)
    
    # Calculate center of the box
    center = np.mean(coords, axis=0)
    
    # Define box parameters
    box_dict = {
        "center_x": center[0],
        "center_y": center[1],
        "center_z": center[2],
        "size_x": size,
        "size_y": size,
        "size_z": size
    }
    
    return box_dict

def clean_pdb_for_gnina(pdb_file, output_folder):
    """
    Writes a chain A receptor PDB for gnina with heteroatoms removed.

    :param pdb_file: Path to input PDB file.
    :param output_folder: Directory where the cleaned PDB should be written.
    :return: Path to the cleaned PDB file.
    """
    pdb_path = Path(pdb_file)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    cleaned_path = output_folder / f"{pdb_path.stem}_chain_A_cleaned.pdb"

    kept_atom_serials = set()
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("ATOM") and line[21] == "A":
                try:
                    kept_atom_serials.add(int(line[6:11]))
                except ValueError:
                    pass

    with open(pdb_path, "r") as in_f, open(cleaned_path, "w") as out_f:
        for line in in_f:
            record = line[:6].strip()

            if record in {"HETATM", "CONECT"}:
                continue

            if record in {"ATOM", "ANISOU", "TER"}:
                if len(line) <= 21 or line[21] != "A":
                    continue

            if record == "ANISOU":
                try:
                    atom_serial = int(line[6:11])
                except ValueError:
                    continue
                if atom_serial not in kept_atom_serials:
                    continue

            out_f.write(line)

    return cleaned_path

def gnina_minimize_defined_box(pdb_file,
                           ligand,
                           output_folder,
                           gnina_path,
                           residues,
                           cnn=None,
                           exhaustiveness=None, 
                           size=None):
    """
    This function both creates poses for the ligand for ESM predictions as well as computes the gnina scores
    for each structure in the pdb_folder. Returns a dataframe with GNina metrics
    
    :param pdb_file: Path to PDB file
    :param ligand: Path to ligand file
    :param output_folder: Output directory for results
    :param gnina_path: Path to gnina executable
    :param residues: List of pocket residues
    :param cnn: CNN model name
    :param exhaustiveness: Exhaustiveness parameter
    :param size: Box size in Angstroms
    """
    if cnn is None:
        cnn = GNINA_CNN_MODEL
    if exhaustiveness is None:
        exhaustiveness = GNINA_EXHAUSTIVENESS
    if size is None:
        size = GNINA_BOX_SIZE
    
    # Create output directory if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Create box dictionary
    box_dict = gnina_box_generator(pdb_file, residues, size)

    # Prepare gnina command
    # Note: Added check to ensure pdb_file is a Path object or string
    pdb_path = Path(pdb_file)
    output_path = Path(output_folder) / f"{pdb_path.stem}_ligand.pdb"
    gnina_command = (
        f'{gnina_path} -r {pdb_file} -l {ligand}  --center_x {box_dict["center_x"]} --center_y {box_dict["center_y"]} --center_z {box_dict["center_z"]} '
        f'--size_x {box_dict["size_x"]} --size_y {box_dict["size_y"]} --size_z {box_dict["size_z"]} '
        f'-o {output_path} '
        f'--exhaustiveness {exhaustiveness} --cnn {cnn}'
    )
    
    # Run the command
    gnina_result = subprocess.run(gnina_command, shell=True, capture_output=True, text=True)
    print(gnina_result)
    # Optional: Debug prints (can comment out for production)
    # print(gnina_command)
    # print(gnina_result.stdout)
    
    output = gnina_result.stdout

    # --- Parsing Logic ---

    # 1. Capture CNN Score
    # Matches: "CNNscore: 0.41463"
    match_cnn_score = re.search(r"CNNscore:\s*([-\d.]+)", output)
    CNN_score = float(match_cnn_score.group(1)) if match_cnn_score else None

    # 2. Capture CNN Affinity
    # Matches: "CNNaffinity: 4.86840"
    match_cnn_aff = re.search(r"CNNaffinity:\s*([-\d.]+)", output)
    CNN_affinity = float(match_cnn_aff.group(1)) if match_cnn_aff else None

    # 3. Capture both Affinity Scores
    # Matches: "Affinity: -6.19290  -0.54008"
    # This regex looks for: "Affinity:", spaces, Group 1 (number), spaces, Group 2 (number)
    match_aff = re.search(r"Affinity:\s*([-\d.]+)\s+([-\d.]+)", output)
    
    if match_aff:
        vina_affinity = float(match_aff.group(1))
        vina_affinity_2 = float(match_aff.group(2))
    else:
        vina_affinity = None
        vina_affinity_2 = None

    return CNN_score, CNN_affinity, vina_affinity, vina_affinity_2

def sample_conformers(molecule, n_conformers, rmsd_cutoff, output, ligand_name="lig"):
    output = Path(output)
    out_dir = output / "conformers"
    
    # Check if conformers already exist
    if out_dir.exists():
        print("Conformers already found. Skipping conformer calculation.")
        # Try to find the lowest energy one based on previous run logic or just return the first found
        # For safety, we assume the best one is formatted as expected if it exists
        stats_file = out_dir / "conformer_stats.csv"
        if stats_file.exists():
             df = pd.read_csv(stats_file)
             # Assuming sorted by energy or index, picking best
             best_idx = df.loc[df['energy'].idxmin()]['conf_id']
             return out_dir / f"{ligand_name}_conf_{int(best_idx)}.pdb"
        return out_dir / f"{ligand_name}_conf_0.pdb" # Fallback

    # Prepare output folder
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load molecule
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {molecule}")
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
### MAIN
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Boltz on a given sequence and ligand.")
    parser.add_argument("--fasta_path", type=str, required=True, help="Path to the input fasta file.")
    parser.add_argument("--smiles_string", type=str, required=True, help="SMILES string of the ligand.")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory to save the Boltz results.")
    parser.add_argument("--ligand_name", type=str, required=True, help="Name of the ligand for output files.")
    parser.add_argument("--num_iterations", type=int, default=1, help="Number of iterations to run Boltz.")
    parser.add_argument("--run_gnina", action="store_true", help="Flag to run Gnina.")
    parser.add_argument("--gnina_pdb_path", type=str, help="Optional receptor PDB path to clean and use for GNINA instead of Boltz-generated PDBs.")
    parser.add_argument("--pocket_residues", nargs="+", type=int, help="List of pocket residues as integers.")
    parser.add_argument("--run_boltz", action="store_true", help="Flag to run Gnina.")
    args = parser.parse_args()

    # Parse fasta file and prepare entries
    fasta_entries = parse_fasta(args.fasta_path)
    
    # Run Boltz for each entry
    for entry in fasta_entries:
        if args.run_boltz:
            # Add file_ID from header for compatibility with yaml generator
            entry['file_ID'] = entry['header']
            entry_output_dir = f"{args.output_dir}/{entry['header']}_{args.ligand_name}"
            os.makedirs(entry_output_dir, exist_ok=True)
            
            for i in range(args.num_iterations):
                iteration_output_dir = f"{entry_output_dir}/iteration_{i}"
                os.makedirs(iteration_output_dir, exist_ok=True)
                
                print(f"Running Boltz iteration {i+1}/{args.num_iterations} for {entry['header']}...")
                result = run_Boltz2(
                    entry, 
                    args.smiles_string, 
                    iteration_output_dir,
                    pocket_list=None,
                    max_dist=BOLTZ_POCKET_MAX_DIST,
                    devices=BOLTZ_DEVICES,
                    sampling_steps=BOLTZ_SAMPLING_STEPS,
                    recycling_steps=BOLTZ_RECYCLING_STEPS,
                    use_msa=BOLTZ_USE_MSA
                )
                
                if result.returncode != 0:
                    print(f"Warning: Boltz failed for {entry['header']} iteration {i+1}")
                else:
                    print(f"Boltz completed successfully for {entry['header']} iteration {i+1}")
            
            # Process Boltz output
            print(f"Processing Boltz output for {entry['header']}...")
            metrics = process_boltz_output(entry_output_dir)
            print(f"Boltz metrics: {metrics}")
        
        # Run GNINA if requested
        if args.run_gnina:
            print(f"Running GNINA for {entry['header']}...")
            entry_output_dir = f"{args.output_dir}/{entry['header']}_{args.ligand_name}"
            gnina_output_dir = f"{entry_output_dir}/gnina_output"
            os.makedirs(gnina_output_dir, exist_ok=True)
            
            # Check if pocket residues are provided
            if not args.pocket_residues:
                print(f"Skipping GNINA for {entry['header']}: no pocket residues provided.")
            else:
                # Generate ligand conformers once (will be cached if already exists)
                print(f"Generating ligand conformers for {entry['header']}...")
                best_conformer_path = sample_conformers(
                    molecule=args.smiles_string, 
                    n_conformers=N_CONFORMERS, 
                    rmsd_cutoff=RMSD_CUTOFF, 
                    output=gnina_output_dir, 
                    ligand_name=LIGAND_NAME
                )
                print(f"Using conformer: {best_conformer_path}")
                
                if args.gnina_pdb_path:
                    cleaned_pdb_path = clean_pdb_for_gnina(args.gnina_pdb_path, gnina_output_dir)
                    pdb_files = [cleaned_pdb_path]
                    print(f"Using cleaned GNINA receptor PDB: {cleaned_pdb_path}")
                else:
                    # Find PDB files generated by Boltz (only from iteration folders)
                    pdb_files = []
                    for i in range(args.num_iterations):
                        iteration_dir = Path(entry_output_dir) / f"iteration_{i}"
                        if iteration_dir.exists():
                            pdb_files.extend(sorted(iteration_dir.rglob('*.pdb')))
                    print(f"Found {len(pdb_files)} PDB files from Boltz iterations")
                
                # Run GNINA on each PDB file
                for pdb_file in pdb_files:
                    try:
                        print(f"Running GNINA on {pdb_file.name}...")
                        cnn_score, cnn_aff, vina_aff, vina_aff_2 = gnina_minimize_defined_box(
                            str(pdb_file),
                            str(best_conformer_path),
                            gnina_output_dir,
                            GNINA_PATH,
                            residues=args.pocket_residues if args.pocket_residues else POCKET_RESIDUES,
                            exhaustiveness=GNINA_EXHAUSTIVENESS,
                            size=GNINA_BOX_SIZE
                        )
                        print(f"GNINA results for {entry['header']}_{args.ligand_name} - CNN Score: {cnn_score}, CNN Affinity: {cnn_aff}, Vina Affinity: {vina_aff}, Vina Affinity 2: {vina_aff_2}")
                    except Exception as e:
                        print(f"Error running GNINA on {pdb_file.name}: {e}")
    
    print("Pipeline completed successfully!")
