"""
Functions for running the script binding_site_designer.py
"""

import warnings
import os
import sys
import json
import ast
import yaml
from Bio import PDB
from Bio.PDB import PDBParser, MMCIFParser, DSSP, NeighborSearch, Selection, Superimposer
from Bio.PDB.Polypeptide import is_aa
from Bio.SVDSuperimposer import SVDSuperimposer
from Bio.SeqUtils import seq1
import numpy as np
import biotite.structure as struc
import biotite.structure.io as bsio
from biotite.structure.io import load_structure as biotite_load_structure

# Biopython loader for generate_redesign_string function

# ============================================================================
# ACTIVE FUNCTIONS FOR BSD PIPELINE (21 functions used by bsd_def.smk)
# ============================================================================

def generate_redesign_string(deNovo_path, reference_path, region=[0,54], chain_id='A', 
                            search_radius=3.0, rmsd_threshold=1.0):
    """
    Identifies residues in de novo structure that differ from wild-type.
    
    Args:
        deNovo_path: Path to RF Diffusion designed structure
        reference_path: Path to original wild-type structure
        region: [start, end] residue numbers for alignment anchor region
        chain_id: Chain ID to analyze (default 'A')
        search_radius: Search radius in Angstroms for finding matching residues (default 3.0)
        rmsd_threshold: RMSD threshold in Angstroms for considering residues as matching (default 1.0)
    
    Returns:
        String of space-separated redesign residues (e.g., "A100 A105 A110")
    """

    WT_structure = load_structure(reference_path)
    dNT_structure = load_structure(deNovo_path)
    
    print("\n" + "="*80)
    print("DEBUG: generate_redesign_string")
    print("="*80)
    print(f"DeNovo path: {deNovo_path}")
    print(f"Reference path: {reference_path}")
    print(f"Search radius: {search_radius} Å")
    print(f"RMSD threshold: {rmsd_threshold} Å")
    print(f"Alignment region: {region[0]} - {region[1]}")
    
    # Aligning Wild Type and de Novo Type by their DNA-binding domain.
    WT_atoms = []
    dNT_atoms = []
    res_index = 1
    aligned_residues = []
    
    for r1, r2 in zip(WT_structure[0]['A'].get_residues(), dNT_structure[0]['A'].get_residues()):
        if r1.get_id()[1] >= region[0] and r1.get_id()[1] < region[1]:
            if is_aa(r1, standard=True) and is_aa(r2, standard=True):
                WT_atoms.append(r1['N'])
                WT_atoms.append(r1['CA'])
                WT_atoms.append(r1['C'])
                WT_atoms.append(r1['O'])
                dNT_atoms.append(r2['N'])
                dNT_atoms.append(r2['CA'])
                dNT_atoms.append(r2['C'])
                dNT_atoms.append(r2['O'])
                aligned_residues.append(r1.get_id()[1])
                res_index += 1
        if res_index >= region[1]:
            break
    
    print(f"\n✓ Aligned {len(aligned_residues)} residues for anchor: {aligned_residues[:5]}...{aligned_residues[-5:]}")
    
    superimpose = Superimposer()
    superimpose.set_atoms(WT_atoms, dNT_atoms)
    alignment_rmsd = superimpose.rms
    print(f"✓ Alignment RMSD: {alignment_rmsd:.3f} Å")
    superimpose.apply(dNT_structure.get_atoms())
    
    # Initiate NeighborSearch object
    WT_atom_space = [atom for atom in WT_structure[0]['A'].get_atoms()]
    NSearch = NeighborSearch(WT_atom_space)
    
    print(f"\n✓ NeighborSearch initialized with {len(WT_atom_space)} WT atoms")
    
    # Retrieve all de Novo Type residues that do not match the Wild Type
    modified_res = ""
    modified_list = []
    matched_list = []
    
    denovo_residues = [res for res in dNT_structure[0]['A'].get_residues() if is_aa(res, standard=True)]
    print(f"\n--- Analyzing {len(denovo_residues)} de novo residues ---\n")
    
    for i, res in enumerate(denovo_residues):
        res_num = res.get_id()[1]
        res_name = res.get_resname()
        
        try:
            proximity_residues = NSearch.search(res['CA'].get_coord(), search_radius, level='R')
        except KeyError:
            print(f"Res {res_num}: ERROR - missing CA atom")
            continue
        
        found_match = False
        
        try:
            res_bb = np.array([res['N'].get_coord(), res['CA'].get_coord(), res['C'].get_coord(), res['O'].get_coord()])
            
            if len(proximity_residues) == 0:
                print(f"Res {res_num} ({res_name}): NO proximity residues found (radius={search_radius}Å)")
            else:
                print(f"Res {res_num} ({res_name}): Found {len(proximity_residues)} proximity residues", end="")
                
                for prox_res in proximity_residues:
                    prox_res_bb = np.array([prox_res['N'].get_coord(), prox_res['CA'].get_coord(), prox_res['C'].get_coord(), prox_res['O'].get_coord()])
                    sup_for_rms = SVDSuperimposer()
                    sup_for_rms.set(res_bb, prox_res_bb)
                    RMS = sup_for_rms.get_init_rms()
                    
                    prox_res_num = prox_res.get_id()[1]
                    prox_res_name = prox_res.get_resname()
                    
                    # Check both RMSD and amino acid type
                    aa_match = res_name == prox_res_name
                    rmsd_match = RMS < rmsd_threshold
                    
                    if rmsd_match and aa_match:
                        print(f" → MATCH with Res {prox_res_num} ({prox_res_name}, RMSD={RMS:.3f}✓)")
                        found_match = True
                        matched_list.append(res_num)
                        break
                    else:
                        details = []
                        if not aa_match:
                            details.append(f"{res_name}→{prox_res_name}")
                        if not rmsd_match:
                            details.append(f"RMSD={RMS:.3f}✗")
                        # print(f"   │ Res {prox_res_num}: {', '.join(details)}", end="")
                
                if not found_match:
                    print(f" → NO GOOD MATCH (checking {len(proximity_residues)} residues)")
        
        except KeyError as key:
            print(f"Res {res_num}: ERROR - missing atom ({key})")
        
        if not found_match:
            modified_res += f" {chain_id}{res_num}"
            modified_list.append(res_num)
            print(f"  ➜ MARKED FOR REDESIGN")
    
    print("\n" + "="*80)
    print(f"SUMMARY:")
    print(f"  - Total residues analyzed: {len(denovo_residues)}")
    print(f"  - Matched (unchanged): {len(matched_list)} residues")
    print(f"  - Modified (for redesign): {len(modified_list)} residues")
    print(f"  - Percent modified: {100*len(modified_list)/len(denovo_residues):.1f}%")
    if len(modified_list) > 0:
        print(f"  - Modified residue numbers: {modified_list[:10]}{'...' if len(modified_list) > 10 else ''}")
    print("="*80 + "\n")
    
    return modified_res.strip()



def process_MPNN_folder(folder,top_n):

    folder = Path(folder)

    sequence_df = pd.DataFrame(columns=["file_ID", "seq", "overall_MPNN_score", "ligand_MPNN_score", "avg_MPNN_score"]) 
    for file in folder.iterdir():
        file_df = process_MPNN_file(file,top_n)
        sequence_df = pd.concat([sequence_df, file_df], ignore_index = True)

    return sequence_df


def process_MPNN_file(file, top_n, select_by='avg_MPNN_score'):
    """Process a single MPNN output file.

    Skips the first header+sequence entry (original sequence) and computes
    an average MPNN score combining overall and ligand confidences.
    """
    # Read headers and sequences
    headers = []
    sequences = []
    with open(file, "r") as faa_file:
        for line in faa_file:
            if line.startswith(">"):
                headers.append(line.strip())
            else:
                sequences.append(line.strip())

    # Parse scores from headers
    overall_scores = []
    ligand_scores = []
    for h in headers:
        m_overall = re.search(r"overall_confidence=([\d.]+)", h)
        m_lig = re.search(r"ligand_confidence=([\d.]+)", h)
        overall_scores.append(float(m_overall.group(1)) if m_overall else np.nan)
        ligand_scores.append(float(m_lig.group(1)) if m_lig else np.nan)

    # If the file contains at least one original entry, drop the first header/sequence
    if len(sequences) > 0:
        # remove first entry (original sequence) from lists if present
        headers = headers[1:]
        sequences = sequences[1:]
        overall_scores = overall_scores[1:]
        ligand_scores = ligand_scores[1:]

    # Store as a dataframe
    sequence_df = pd.DataFrame(columns=["file_ID", "seq", "overall_MPNN_score", "ligand_MPNN_score", "avg_MPNN_score"]) 
    file_name = str(file.name)
    # Remove common extension if present (e.g., .fa)
    if file_name.lower().endswith('.fa'):
        file_name = file_name[:-3]
    elif file_name.lower().endswith('.fasta'):
        file_name = file_name[:-6]

    for idx, seq in enumerate(sequences, start=1):
        entry_name = f"{file_name}_seq_{idx}"
        overall = overall_scores[idx-1] if idx-1 < len(overall_scores) else np.nan
        ligand = ligand_scores[idx-1] if idx-1 < len(ligand_scores) else np.nan
        # Compute average (ignore NaNs when possible)
        avg = float(np.nanmean([overall, ligand])) if (not np.isnan(overall) or not np.isnan(ligand)) else np.nan
        row = pd.Series([entry_name, seq, overall, ligand, avg], index=["file_ID", "seq", "overall_MPNN_score", "ligand_MPNN_score", "avg_MPNN_score"])
        sequence_df.loc[len(sequence_df)] = row

    # Sort by desired metric (descending) and keep top N
    if select_by not in sequence_df.columns:
        select_by = 'avg_MPNN_score'
    top_df = sequence_df.sort_values(by=select_by, ascending=False).head(top_n).reset_index(drop=True)

    return top_df

### Folding models ###############################################################################################################################


def run_AlphaFold3(row, output_path,msa_flag=False,ligand_SMILES=None, env=None):
    # Create temporary JSON file
    if msa_flag:
        AF_json_generator_w_MSA(row, output_path)
        AF_json_generator_cofold_w_MSA(row,ligand_SMILES, output_path)
    else:
        AF_json_generator_wo_MSA(row, output_path)
        AF_json_generator_cofold(row,ligand_SMILES, output_path)
    # Define the command
    af_command = (
        f"conda run -p /home/ingemar/anaconda3/envs/alphafold3 python /mnt/data/alphafold3/run_alphafold.py "
        f"--json_path={output_path}/temp.json "
        f"--output_dir={output_path}/AF3_prediction "
        f"--db_dir=/mnt/data/alphafold3/alphafold_databases/ "
        f"--model_dir=/mnt/data/alphafold3/models/ "
        f"--num_diffusion_samples=1"
    )

    af_command_cofold = (
        f"conda run -p /home/ingemar/anaconda3/envs/alphafold3 python /mnt/data/alphafold3/run_alphafold.py "
        f"--json_path={output_path}/temp_cofold.json "
        f"--output_dir={output_path}/AF3_cofold "
        f"--db_dir=/mnt/data/alphafold3/alphafold_databases/ "
        f"--model_dir=/mnt/data/alphafold3/models/ "
        f"--num_diffusion_samples=1"
    )

    # Clone the current environment and add the JAX flag
    current_env = os.environ.copy()
    current_env["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"

    # Execute
    if env:
        current_env.update(env)

    subprocess.run(af_command, shell=True, env=current_env)
    subprocess.run(af_command_cofold, shell=True, env=current_env)


def run_chai(row, output_path, env=None):
    # Ensure output path exists
    Path(output_path).mkdir(parents=True, exist_ok=True)
    # 
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
    # Run chai
    chai_command =f" conda run -n chai chai-lab fold {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True, env=env)

    return 


def run_chai_w_MSA(row, output_path, env=None):
    #
    Path(output_path).mkdir(parents=True, exist_ok=True)
    #
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
    # Run chai
    chai_command =f" conda run -n chai chai-lab fold --use-msa-server --use-templates-server {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True, env=env)

    return 


def run_chai_cofold(row, output_path,ligand_smiles, env=None):
    #
    Path(f"{output_path}").mkdir(parents=True, exist_ok=True)
    #
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
        f.write(f">ligand|ligand\n")      # Header line
        f.write(f"{ligand_smiles}\n")         # The actual ligand SMILES
    # Run chai
    chai_command =f" conda run -n chai chai-lab fold {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True, env=env)

    return 


def run_chai_cofold_w_MSA(row, output_path,ligand_smiles, env=None):
    #
    Path(f"{output_path}").mkdir(parents=True, exist_ok=True)
    #
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
        f.write(f">ligand|ligand\n")      # Header line
        f.write(f"{ligand_smiles}\n")         # The actual ligand SMILES
    # Run chai
    chai_command =f" conda run -n chai chai-lab fold --use-msa-server --use-templates-server {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True, env=env)

    return 


def boltz_yaml_generator(row, yaml_path, ligand_smiles, pocket_list, max_dist=5.0):
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
                        "sequence": seq,
                        "msa": "empty"  
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
                        "sequence": seq,
                        "msa": "empty"  
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


def get_centroid_distance(structure_path, ligand_id, site_residues, center_cutoff):
    """
    Helper function: Parses PDB and calculates centroid distance.
    Returns: (is_valid, distance)
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('struct', structure_path)
    except Exception as e:
        print(f"Error parsing {structure_path}: {e}")
        return False, np.nan

    # 1. Extract Ligand Atoms
    ligand_atoms = []
    for res in Selection.unfold_entities(structure, 'R'):
        # Check against ligand_id (adjust if your PDB uses different naming)
        if res.get_resname() == ligand_id:
            ligand_atoms.extend(res.get_atoms())
            
    if not ligand_atoms:
        return False, np.nan  # Ligand missing

    # 2. Extract Active Site Atoms
    site_atoms = []
    for res in Selection.unfold_entities(structure, 'R'):
        if res.id[1] in site_residues:
            site_atoms.extend(res.get_atoms())

    if not site_atoms:
        return False, np.nan  # Active site missing

    # 3. Calculate Centroids using NumPy
    lig_coords = np.array([atom.get_coord() for atom in ligand_atoms])
    site_coords = np.array([atom.get_coord() for atom in site_atoms])
    
    lig_centroid = np.mean(lig_coords, axis=0)
    site_centroid = np.mean(site_coords, axis=0)

    # 4. Calculate Distance
    distance = np.linalg.norm(lig_centroid - site_centroid)
    is_valid = distance <= center_cutoff
    
    return is_valid, distance


def check_cofold_validity(df, path_to_structures, ligand_id, site_residues, center_cutoff=4.0,extension=".pdb"):
    """
    Updates the dataframe with 'distance' and 'is_valid' columns based on ligand placement.
    """
    distances = []
    validity = []
    
    print(f"Processing {len(df)} structures for ligand '{ligand_id}'...")

    for index, row in df.iterrows():
        file_id = row['file_ID'].split('.')[0]  # Remove extension if present
        file_id = file_id.lower()  # Ensure lowercase for consistency
        # Construct full path (handles missing extension if needed)
        structure_path = os.path.join(path_to_structures, f"{file_id}{extension}")
        
        if not os.path.exists(structure_path):
            print(f"File not found: {structure_path}")
            distances.append(np.nan)
            validity.append(False)
            continue
            
        # Call the helper
        is_valid, dist = get_centroid_distance(structure_path, ligand_id, site_residues, center_cutoff)
        
        distances.append(dist)
        validity.append(is_valid)

    # Assign new columns
    df['distance'] = distances
    df['is_valid'] = validity
    
    return df
### 3D FILTERING #######################################################################################################################################


def filter_dataframe_1D(df, window_size, threshold, seq_col='seq', aa= 'A'):
    """
    Removes rows from the DataFrame where the sequence contains 
    too many Alanines in any sliding window.

    Args:
        df (pd.DataFrame): Input dataframe.
        window_size (int): Size of the window (y).
        threshold (int): Max allowed Alanines in a window (x).
        seq_col (str): The name of the sequence column. Defaults to 'seq'.
        
    Returns:
        pd.DataFrame: A filtered copy of the dataframe.
    """
    # Create a mask: True if the sequence HAS the bad region
    # We use a lambda to pass the extra arguments to your function
    mask = df[seq_col].apply(lambda s: sliding_window_1D_filter(s, window_size, threshold,aa))
    
    # We keep the rows where the mask is False (using the bitwise NOT operator ~)
    df_filtered = df[~mask].copy()
    
    return df_filtered


def threed_params_1_df(folder, output_folder,output_name, original_path, gnina_pocket_residues, clash_distance=2.0, bond_distance=1.2,
                       ligand_path=None, gnina_path="gnina", cnn="default",
                       exhaustiveness=8, gnina_box_size=4):
    """
    Given a folder with structures, returns the following df: "file_ID", "sequence", "pLDDT mean", "pLDDT std dev", "RMSD", "TMscore", "clashes","Gnina_CNNscore","Gnina_affinity"
    """
    # Skip if output already exists and is not empty
    if os.path.exists(f"{output_folder}/{output_name}") and os.path.getsize(f"{output_folder}/{output_name}") > 0:
        print(f"3D metrics already exists at {output_folder}/{output_name}. Skipping generation.")
        threed_df = pd.read_csv(f"{output_folder}/{output_name}")
        return threed_df
    # Initialise df
    elements = ["file_ID", "sequence", "pLDDT_mean", "pLDDT_std", "RMSD", "TMscore", "num_clashes","clashes_per_atom","Gnina_CNNscore","Gnina_affinity"]
    threed_df = pd.DataFrame(columns = elements)

    # Parse folder
    folder = Path(folder)
    for file in folder.iterdir():
        file_path = str(file)
        if file_path.endswith(".pdb") or file_path.endswith(".cif"):
            # pLDDT
            pLDDT_mean,pLDDT_stdev= extract_plddt(file_path)
            # RMSD
            RMSD, TMscore = batch_tm_align(original_path,file_path)
            # Clashes
            num_clashes,clash_atom = detect_clashes(file_path, clash_distance=clash_distance, bond_distance=bond_distance)
            # Sequence
            seq = extract_sequence(file_path, chain_id='A')["coordinate_sequence"]
            # GNINA scores
            Gnina_CNNscore, Gnina_affinity,vina_affinity,vina_affinity_2 = gnina_minimize_defined_box(pdb_file=file,
                                                                    ligand=ligand_path,
                                                                    output_folder=f"{output_folder}/ESM_Gnina_minimized/",
                                                                    gnina_path=gnina_path,
                                                                    cnn=cnn,
                                                                    exhaustiveness=exhaustiveness,
                                                                    residues=gnina_pocket_residues,
                                                                    size=gnina_box_size)
            # Generate row
            row = pd.Series([file.name, seq, pLDDT_mean, pLDDT_stdev, RMSD, TMscore, num_clashes,clash_atom,Gnina_CNNscore,Gnina_affinity], index=list(threed_df.columns))
            threed_df.loc[len(threed_df)] = row

    print(threed_df.head(5))
    threed_df.to_csv(f"{output_folder}/{output_name}", index=False)
    return threed_df


def global_score(df, weights={"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2}):
    # Iterate through the weights to normalize each relevant column
    for col, weight in weights.items():
        col_min = df[col].min()
        col_max = df[col].max()
        
        # Create a new column with suffix '_norm'
        # Handle the edge case where all values are identical (div by zero)
        if col_max == col_min:
            df[f"{col}_norm"] = 1.0 
        else:
            df[f"{col}_norm"] = (df[col] - col_min) / (col_max - col_min)

    # Calculate the weighted sum using the newly created normalized columns
    df['global_score'] = sum(df[f"{col}_norm"] * weight for col, weight in weights.items())
    
    return df


def threed_filter_1_df(df,output_folder, weights, min_plddt=0.8, min_rmsd=1, max_rmsd=8, max_clashes=1.01, top_n_score=10,top_n_gnina=10):

    # Skip if output already exists and is not empty
    if os.path.exists(f"{output_folder}/ESM_filtered_df.csv") and os.path.getsize(f"{output_folder}/ESM_filtered_df.csv") > 0:
        print(f"3D metrics already exists at {output_folder}/ESM_filtered_df.csv. Skipping generation.")
        df = pd.read_csv(f"{output_folder}/ESM_filtered_df.csv")
        return df
    
    # Drop rows with low pLDDT
    df = df[df["pLDDT_mean"] >= min_plddt]

    # Drop rows with RMSD out of range
    df = df[(df["RMSD"] >= min_rmsd) & (df["RMSD"] <= max_rmsd)]

    # Drop rows with too many clashes
    df = df[df["clashes_per_atom"] <= max_clashes]

    # Global score calculation
    df = global_score(df, weights)

    # Select top n by global score
    sorted_df = df.sort_values(by="global_score", ascending=False).reset_index(drop=True)
    sorted_df = sorted_df.head(top_n_score)

    # Select top n by GNINA
    sorted_df = sorted_df.sort_values(by="Gnina_affinity", ascending=True).reset_index(drop=True)
    sorted_df = sorted_df.head(top_n_gnina)
    sorted_df = sorted_df.sort_values(by="global_score", ascending=False).reset_index(drop=True)

    return sorted_df


def threed_filter_2_df(df_list, output_folder, output_names, weights, min_plddt=0.8, min_rmsd=1, max_rmsd=8, max_clashes=1.01, top_n_score=10, top_n_gnina=10):
    filtered_dfs = []
    
    # 2. First Pass: Apply filters and calculate individual global scores
    for df in df_list:
        # Clean file ID
        df["file_ID"] = df["file_ID"].str.split('.').str[0]
        
        # Filter by physical constraints
        df = df[
            (df["pLDDT_mean"] >= min_plddt) & 
            (df["RMSD"] >= min_rmsd) & 
            (df["RMSD"] <= max_rmsd) & 
            (df["clashes_per_atom"] <= max_clashes)
        ].copy()
        
        # Calculate global score for this specific replicate
        df = global_score(df, weights)
        filtered_dfs.append(df)

    # 3. Find intersection of IDs (only keep IDs that survived filters in ALL dataframes)
    if not filtered_dfs:
        return []
        
    common_ids = set(filtered_dfs[0]["file_ID"])
    for df in filtered_dfs[1:]:
        common_ids.intersection_update(set(df["file_ID"]))
    
    if not common_ids:
        print("No designs survived the filters across all dataframes.")
        return []

    # 4. Build the intermediate "temp_df" to sum scores
    # Initialize with the common IDs
    temp_df = pd.DataFrame({"file_ID": list(common_ids)})
    temp_df["joint_global_score"] = 0.0
    temp_df["joint_gnina"] = 0.0

    # Add up the scores from each filtered dataframe
    for df in filtered_dfs:
        # Filter this df to only common IDs to ensure alignment
        subset = df[df["file_ID"].isin(common_ids)].set_index("file_ID")
        
        # Map scores to the temp_df
        temp_df["joint_global_score"] += temp_df["file_ID"].map(subset["global_score"]).fillna(0)
        temp_df["joint_gnina"] += temp_df["file_ID"].map(subset["Gnina_affinity"]).fillna(0)

    # 5. Ranking Selection on the Joint Scores
    temp_df.to_csv(f"{output_folder}/final_scores.csv", index=False)
    # Sort by Joint Global Score -> Keep Top N
    temp_df = temp_df.sort_values(by="joint_global_score", ascending=False).head(top_n_score)
    temp_df.to_csv(f"{output_folder}/final_scores_topn_filtered.csv", index=False)
    # Sort by Joint Gnina (ascending usually for affinity) -> Keep Top N
    temp_df = temp_df.sort_values(by="joint_gnina", ascending=True).head(top_n_gnina)
    temp_df.to_csv(f"{output_folder}/final_scores_topgnina_filtered.csv", index=False)
    
    # 6. Final Filter: Return the original dataframes containing only the winning IDs
    winning_ids = temp_df["file_ID"].values
    
    final_output_dfs = []
    for df, name in zip(filtered_dfs, output_names):
        final_df = df[df["file_ID"].isin(winning_ids)].copy()
        
        # Optional: Save individual CSVs for each replicate
        final_df.to_csv(f"{output_folder}/{name}.csv", index=False)
        final_output_dfs.append(final_df)

    return final_output_dfs

### VISUALIZATION #######################################################################################################################################

def plot_four_scatters_with_region(dfs, titles, x_col, y_col, x_range, y_range, xlim=None, ylim=None, output_file=None):
    """
    Plots a 2x2 grid of scatter plots, overlays a rectangle on each, 
    and displays the count of designs within the target region.
    """
    
    # Create a 2x2 grid of subplots
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten() # Flattens the 2x2 array into a 1D list of 4 axes for easy looping
    
    for idx, (df, title) in enumerate(zip(dfs, titles)):
        ax = axes[idx]
        
        # 1. Create the scatter plot
        # Use .empty and check if columns exist to avoid errors if a model failed and returned an empty df
        if not df.empty and x_col in df.columns and y_col in df.columns:
            ax.scatter(df[x_col], df[y_col], alpha=0.6, label='Data Points')
            
            # 4. Count designs in region using the external function
            design_count = count_designs_in_region(df, x_col, y_col, x_range, y_range)
        else:
            design_count = 0
            ax.text(0.5, 0.5, "No Data", ha='center', va='center', transform=ax.transAxes)

        # 2. Calculate rectangle parameters
        width = x_range[1] - x_range[0]
        height = y_range[1] - y_range[0]
        
        # 3. Create the rectangle patch
        rect = patches.Rectangle(
            (x_range[0], y_range[0]), 
            width, height, linewidth=2, 
            edgecolor='red', facecolor='red', alpha=0.2, label='Target Region'
        )
        ax.add_patch(rect)
        
        # Add text box with the count
        text_str = f'Designs in region: {design_count}'
        props = dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray')
        ax.text(0.05, 0.95, text_str, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=props)
        
        # 5. Apply Axis Limits
        if xlim is not None: ax.set_xlim(xlim)
        if ylim is not None: ax.set_ylim(ylim)
        
        # Formatting
        ax.set_xlabel(x_col)
        ax.set_ylabel(y_col)
        ax.set_title(title)
        ax.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout() # Prevents overlapping of titles and labels
    
    # Save or Show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.show()
        
    plt.close(fig)


def plot_comparative_distributions(df_list, columns, df_labels=None, cols=2, output_file=None):
    """
    Plots distribution curves for specific columns across multiple dataframes.
    
    df_list: List of pandas DataFrames
    columns: List of column names to plot
    df_labels: List of names for the dataframes (for the legend)
    cols: Number of columns in the subplot grid
    output_file: Path to save the image for Snakemake
    """
    import math
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    # 1. Setup Labels
    if df_labels is None:
        df_labels = [f'Dataset {i+1}' for i in range(len(df_list))]
        
    if len(df_labels) != len(df_list):
        raise ValueError("Length of df_labels must match length of df_list")

    # 2. Setup Grid
    n_plots = len(columns)
    
    # --- FIX 1: Don't force 2 columns if we only have 1 plot! ---
    actual_cols = min(cols, n_plots) 
    rows = math.ceil(n_plots / actual_cols)
    
    fig, axes = plt.subplots(rows, actual_cols, figsize=(actual_cols * 6, rows * 4), constrained_layout=True)
    
    # --- FIX 2: Safely flatten regardless of grid size ---
    if hasattr(axes, 'flatten'):
        axes_flat = axes.flatten()
    else:
        axes_flat = [axes]

    # 3. Plotting Loop
    for i, col_name in enumerate(columns):
        ax = axes_flat[i]
        
        # Iterate through every dataframe and add its curve to the current plot
        for j, df in enumerate(df_list):
            if col_name in df.columns:
                sns.kdeplot(
                    data=df, 
                    x=col_name, 
                    label=df_labels[j], 
                    ax=ax, 
                    fill=True, 
                    alpha=0.1
                )
        
        ax.set_title(f'Distribution of {col_name}')
        ax.set_xlabel(col_name)
        ax.set_ylabel('Density')
        ax.legend(title="Data Source")
        ax.grid(True, linestyle=':', alpha=0.5)

    # 4. Cleanup
    for k in range(i + 1, len(axes_flat)):
        axes_flat[k].axis('off')
        
    # --- NEW: Save logic for Snakemake ---
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
    else:
        plt.show()
        
    plt.close(fig)
### AUXILIARY FUNCTIONS ###############################################################################################################################

def process_chai_folder(base_path, output_base_path):
    """
    Scans for 'CHAI_prediction' and 'CHAI_cofold' folders.
    Inside each sequence subfolder, it finds 'pred.model_idx_0.cif',
    renames it to the sequence ID (folder name), and moves it to 
    'CHAI_prediction_pdbs' or 'CHAI_cofold_pdbs'.

    Args:
        base_path (str): The root path containing CHAI_prediction/CHAI_cofold folders.
        output_base_path (str): The destination root path for the new _pdbs folders.
    """
    base_path = Path(base_path)
    output_base_path = Path(output_base_path)

    # Map source folders to destination folders
    folders_to_process = {
        "CHAI_prediction": "CHAI_prediction_pdbs",
        "CHAI_cofold": "CHAI_cofold_pdbs"
    }

    for source_name, dest_name in folders_to_process.items():
        source_dir = base_path / source_name
        dest_dir = output_base_path / dest_name

        if not source_dir.exists():
            continue

        print(f"Processing {source_name} -> {dest_name}...")
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Iterate over the subdirectories (each represents a sequence/structure)
        # Example subdir: structure_pos0_conf0_rfdesign_1_seq_1.pdb
        for subdir in source_dir.iterdir():
            if subdir.is_dir():
                # Define the specific file we want (Model 0)
                model_0_file = subdir / subdir.name.lower() / "pred.model_idx_0.cif"
                print (model_0_file)

                if model_0_file.exists():
                    try:
                        # The renaming logic:
                        # Use the folder name as the ID. 
                        # If the folder name ends in '.pdb' (as seen in the image), we strip it to be clean.
                        folder_name = subdir.name
                        if folder_name.endswith('.pdb'):
                            new_name = folder_name[:-4] + ".cif"
                        else:
                            new_name = folder_name + ".cif"

                        dest_file = dest_dir / new_name
                        print(f"  Moving {model_0_file} to {dest_file}...")

                        # Copy and rename
                        shutil.copy2(model_0_file, dest_file)
                        
                    except Exception as e:
                        print(f"  Error moving {model_0_file}: {e}")
                else:
                    # Optional: Warn if model 0 is missing in a folder
                    # print(f"  Warning: pred.model_idx_0.cif not found in {subdir.name}")
                    pass

    return


def process_AF3_folder(base_path, output_base_path):
    """
    Scans the base directory for 'AF3_prediction' and 'AF3_cofold' folders.
    Moves/Copies the resulting .cif files from the subdirectories into 
    flat 'AF3_prediction_pdbs' and 'AF3_cofold_pdbs' folders.
    
    Args:
        base_path (str): The root path containing AF3_prediction/AF3_cofold folders.
        output_base_path (str): The destination root path for the new _pdbs folders.
    """
    base_path = Path(base_path)
    output_base_path = Path(output_base_path)
    
    # Define the mapping of source folder names to destination folder names
    # Key: source subfolder name, Value: destination subfolder name
    folders_to_process = {
        "AF3_prediction": "AF3_prediction_pdbs",
        "AF3_cofold": "AF3_cofold_pdbs"
    }

    for source_name, dest_name in folders_to_process.items():
        source_dir = base_path / source_name
        dest_dir = output_base_path / dest_name
        
        # Skip if the source directory doesn't exist (e.g. if you didn't run cofold)
        if not source_dir.exists():
            continue

        print(f"Processing {source_name} -> {dest_name}...")
        
        # Create destination directory
        dest_dir.mkdir(parents=True, exist_ok=True)
        
        # AF3 outputs are usually nested: Source_Dir / Job_Name / model.cif
        # We search recursively for .cif files
        cif_files = list(source_dir.rglob("*.cif"))
        
        if not cif_files:
            print(f"  No .cif files found in {source_dir}")
            continue
            
        count = 0
        for cif_file in cif_files:
            try:
                # The file is usually named 'model.cif' or similar generic names inside a named folder.
                # To avoid collisions and keep the sequence ID, we prefer the name of the PARENT folder 
                # (which is usually the job name/sequence header).
                
                # Example: .../structure_pos0_..._seq_1.pdb/structure_..._model.cif
                # We use the filename itself if it's descriptive, otherwise fall back to parent.
                
                new_filename = cif_file.name
                
                # If the filename is just "model.cif", prepend the parent folder name
                if new_filename == "model.cif":
                    new_filename = f"{cif_file.parent.name}.cif"
                
                dest_file = dest_dir / new_filename
                
                # Copy the file (using copy2 to preserve metadata)
                shutil.copy2(cif_file, dest_file)
                count += 1
                
            except Exception as e:
                print(f"  Error moving {cif_file}: {e}")
                
        print(f"  Successfully moved {count} files to {dest_dir}")

    return




def process_Boltz_folder(input_folder, output_pdbs_folder):
    """
    Parses a Boltz prediction folder, moves PDBs to a central folder, 
    and aggregates confidence and affinity metrics into a DataFrame.

    Args:
        input_folder (str): Path to the root BOLTZ_prediction folder.
        output_pdbs_folder (str): Path where the resulting PDBs should be moved/copied.

    Returns:
        pd.DataFrame: A dataframe containing file IDs, confidence metrics, and affinity metrics.
    """
    input_path = Path(input_folder)
    output_pdb_path = Path(output_pdbs_folder)
    
    # Create output directory if it doesn't exist
    output_pdb_path.mkdir(parents=True, exist_ok=True)

    data_rows = []
    
    # Recursively find all confidence JSON files
    json_files = list(input_path.rglob("confidence_*.json"))
    
    if not json_files:
        print(f"No confidence JSON files found in {input_folder}")
        return pd.DataFrame()


    for json_file in json_files:
        try:
            # 1. Parse the Confidence JSON file
            with open(json_file, 'r') as f:
                data = json.load(f)

            # Derive file_ID (e.g., structure_..._model_0)
            file_id = json_file.stem.replace("confidence_", "")
            file_id = file_id.replace("_model_0", "")
            # Extract standard confidence fields
            row = {
                "file_ID": file_id,
                "confidence_score": data.get("confidence_score"),
                "ptm": data.get("ptm"),
                "iptm": data.get("iptm"),
                "ligand_iptm": data.get("ligand_iptm"),
                "protein_iptm": data.get("protein_iptm"),
                "complex_plddt": data.get("complex_plddt"),
                "complex_iplddt": data.get("complex_iplddt"),
                "complex_pde": data.get("complex_pde"),
                "complex_ipde": data.get("complex_ipde")
            }
            
            # 2. Look for corresponding Affinity JSON file
            # Assuming naming convention: affinity_[file_ID].json
            affinity_file = json_file.parent / f"affinity_{file_id}.json"
            
            if affinity_file.exists():
                try:
                    with open(affinity_file, 'r') as af:
                        affinity_data = json.load(af)
                        
                        # Add specific affinity metrics to the row
                        # We use .get() to avoid errors if keys are missing
                        row["affinity_pred_value"] = affinity_data.get("affinity_pred_value")
                        row["affinity_probability_binary"] = affinity_data.get("affinity_probability_binary")
                        
                        # Add ensemble member metrics if they exist (e.g., value1, value2)
                        # We iterate the keys to capture dynamic fields like value1, value2, etc.
                        for key, value in affinity_data.items():
                            if key not in row: # Avoid overwriting existing keys if any overlap
                                row[key] = value
                                
                except Exception as e:
                    print(f"Error processing affinity file {affinity_file}: {e}")
            else:
                # If you expect every file to have affinity, you might want to warn here.
                # print(f"Warning: Affinity file not found for {file_id}")
                pass

            # 3. Locate and Move the corresponding PDB file
            expected_pdb_name = f"{file_id}.pdb"
            source_pdb_path = json_file.parent / expected_pdb_name
            
            if source_pdb_path.exists():
                dest_pdb_path = output_pdb_path / expected_pdb_name
                shutil.copy2(source_pdb_path, dest_pdb_path)
            else:
                # Try finding it with the original filename logic just in case
                expected_pdb_name_alt = json_file.name.replace("confidence_", "").replace(".json", ".pdb")
                source_pdb_path_alt = json_file.parent / expected_pdb_name_alt
                
                if source_pdb_path_alt.exists():
                    shutil.copy2(source_pdb_path_alt, output_pdb_path / expected_pdb_name_alt)
                else:
                    print(f"Warning: Corresponding PDB not found for {json_file.name}")

            data_rows.append(row)
            
        except Exception as e:
            print(f"Error processing {json_file}: {e}")

    # Create DataFrame
    df = pd.DataFrame(data_rows)
    
    # Reorder columns to ensure file_ID is first
    if not df.empty and 'file_ID' in df.columns:
        cols = ['file_ID'] + [c for c in df.columns if c != 'file_ID']
        df = df[cols]

    return df


