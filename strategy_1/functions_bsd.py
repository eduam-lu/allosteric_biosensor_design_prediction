"""
Functions for running the script binding_site_designer.py
"""

import warnings
import os
import sys
import json
from Bio import PDB
from Bio.PDB import PDBParser, DSSP, NeighborSearch
from Bio.SeqUtils import seq1

import biotite.structure as struc
import biotite.structure.io as bsio
from biotite.structure.io import load_structure
import numpy as np
from tmtools import tm_align

import pandas as pd
import re
from pymol import cmd 
from pathlib import Path
import subprocess



### PDB INFO EXTRACTION ###########################################################################################################################

def extract_pdb_info(pdb_path, 
                    first_shell_distance=5.0, second_shell_distance=5.0, 
                    user_defined_active_site=None, user_defined_residues=None):
    # Extract chain id
    chain_id = get_chain_id(pdb_path)
    # Extract sequence
    sequence_dict = extract_sequence(pdb_path, chain_id)
    # Extract DSSP
    dssp_string = extract_dssp_string(pdb_path)
    # Define active site
    active_site_dict = define_active_site(pdb_path, 
                                          first_shell_distance=5.0, second_shell_distance=5.0, 
                                          user_defined_active_site=None, user_defined_residues=None)
    # Prepare final dictionary
    pdb_info = {
        "pdb_id": sequence_dict["pdb_id"],
        "chain_id": chain_id,
        "sequence_seqres": sequence_dict["seqres"],
        "sequence_coordinate": sequence_dict["coordinate_sequence"],
        "residue_count": sequence_dict["residue_count"],
        "found_residues": sequence_dict["found_residues"],
        "missing_residues": sequence_dict["missing_residues"],
        "missing_positions": sequence_dict["missing_positions"],
        "start_residue_number": sequence_dict["start_residue_number"],
        "sequence_length": sequence_dict["sequence_length"],
        "dssp_string": dssp_string,
        "active_site": active_site_dict["active_site"],
        "first_shell": active_site_dict["first_shell"],
        "second_shell": active_site_dict["second_shell"],
        "user_defined_active_site": active_site_dict["user_defined_active_site"],
        "user_defined_residues": active_site_dict["user_defined_residues"]
    }
    printable_pdb_info = {
        "pdb_id": sequence_dict["pdb_id"],
        "chain_id": chain_id,
        "sequence_seqres": sequence_dict["seqres"],
        "sequence_coordinate": sequence_dict["coordinate_sequence"],
        "residue_count": str(sequence_dict["residue_count"]),
        "found_residues": str(sequence_dict["found_residues"]),
        "missing_residues": str(sequence_dict["missing_residues"]),
        "missing_positions": str(sequence_dict["missing_positions"]),
        "start_residue_number": str(sequence_dict["start_residue_number"]),
        "sequence_length": str(sequence_dict["sequence_length"]),
        "dssp_string": dssp_string,
        "active_site": str(active_site_dict["active_site"]),
        "first_shell": str(active_site_dict["first_shell"]),
        "second_shell": str(active_site_dict["second_shell"]),
        "user_defined_active_site": str(active_site_dict["user_defined_active_site"]),
        "user_defined_residues": str(active_site_dict["user_defined_residues"])
    }
    return pdb_info, printable_pdb_info

def get_chain_id(pdb_path):
    """
    Parses a PDB file and returns the chain ID.
    If multiple chains exist, it returns the first one found.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return None # Or raise the error depending on your preference

    # Access the first model (Model 0)
    if len(structure) > 0:
        model = structure[0]
        # Get iterator of chains and convert to list to access by index
        chains = list(model.get_chains())
        
        if chains:
            return chains[0].id
            
    return None

def extract_sequence(pdb_path, chain_id='A'):
    """
    Main wrapper to parse file and combine info.
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {"error": str(e)}

    # 1. Get SEQRES sequence
    seqres = get_seqres_sequence(structure, chain_id)
    
    # 2. Get Coordinate Sequence
    if chain_id in structure[0]:
        coord_seq, res_count = get_coordinate_sequence(structure[0][chain_id])
    else:
        coord_seq, res_count = ("", 0)
    
    if seqres is None:
        seqres = coord_seq  # Fallback to coordinate sequence if SEQRES not found
    
    # 3. Complete C term missing residues with '-'
    if len(coord_seq) < len(seqres):
        coord_seq += '-' * (len(seqres) - len(coord_seq))
    
    missing_residues = coord_seq.count('-')
    found_residues = len(coord_seq) - missing_residues
    missing_positions = [i+1 for i, res in enumerate(coord_seq) if res == '-']

    # 4. Find first position that's not - to determine starting residue number
    first_residue_pos = next((i for i, res in enumerate(coord_seq) if res != '-'), None)
    if first_residue_pos is not None:
        start_residue_number = first_residue_pos + 1  # +1 for 1-based indexing
    else:
        start_residue_number = None  # All residues are missing

    return {
        "pdb_id": str(pdb_path).split('/')[-1].split('.')[0],
        "seqres": seqres,
        "coordinate_sequence": coord_seq,
        "residue_count": len(coord_seq),
        "found_residues": found_residues,
        "missing_residues": missing_residues,
        "missing_positions": missing_positions,
        "start_residue_number": start_residue_number,
        "sequence_length": len(coord_seq)
    }

def get_seqres_sequence(structure, chain_id='A'):
    """
    Retrieves the SEQRES sequence from the PDB header.
    Returns None if not found.
    """
    if 'seqres' in structure.header:
        seqres_dict = structure.header['seqres']
        if chain_id in seqres_dict:
            # Biopython stores this as a sequence string
            return str(seqres_dict[chain_id])
    return None

def get_coordinate_sequence(chain):
    """
    Extracts sequence from atom coordinates, accounting for:
    1. Missing N-term residues (pads with '-' if start > 1).
    2. Missing internal residues (pads with '-' for missing numbers).
    
    Returns:
        tuple: (sequence_string, residue_count)
    """
    # Filter for standard amino acids (ignore waters/ligands)
    residues = [res for res in chain if PDB.is_aa(res)]
    
    if not residues:
        return "", 0

    # Extract residue numbers (id[1]) and 1-letter codes
    # res.id looks like (' ', 10, ' ') -> we want the 10
    res_map = {res.id[1]: seq1(res.resname) for res in residues}
    
    start_res_num = min(res_map.keys())
    end_res_num = max(res_map.keys())
    
    # Check if we need to pad the beginning (N-term missing)
    # User requirement: Check if first residue is not 1
    sequence_parts = []
    
    # Iterate from 1 (standard start) to the last observed residue
    # This automatically handles N-term gaps (1 to start) and internal gaps
    for i in range(1, end_res_num + 1):
        if i in res_map:
            sequence_parts.append(res_map[i])
        else:
            sequence_parts.append('-') # Gap detected
            
    final_seq = "".join(sequence_parts)
    residue_count = len(residues) # Count of actual physical residues
    
    return final_seq, residue_count

def extract_dssp_string(input_file):
    # Parse the PDB file
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_file)
    model = structure[0]  # DSSP works on one model, usually model 0

    # DSSP calculation
    # Ensure that you have DSSP installed and correctly configured
    dssp = DSSP(model, input_file)

    # Get DSSP codes (secondary structure assignments)
    dssp_codes = ''.join([dssp[key][2] for key in dssp.keys()])

    return dssp_codes

def define_active_site(pdb_path,user_defined_active_site=None, user_defined_residues = None, first_shell_distance=5.0, second_shell_distance=5.0):
    # If active site already defined, return it. Ensures desired residues remain too
    if user_defined_active_site is not None:
        active_site = user_defined_active_site
        for i in user_defined_residues:
            if i not in active_site:
                active_site.append(i)
    
    # Pymol selection based on distance to ligand
    first_shell = detect_first_shell(pdb_path, distance=first_shell_distance)
    second_shell = detect_second_shell(pdb_path, first_shell, distance=second_shell_distance)
    active_site = first_shell

    # Include user-defined residues
    if user_defined_residues is not None:
        for i in user_defined_residues:
            if i not in active_site:
                active_site.append(i)

    return {"active_site": active_site,
            "first_shell": first_shell,
            "second_shell": second_shell,
            "user_defined_active_site": user_defined_active_site,
            "user_defined_residues": user_defined_residues}

def detect_first_shell(pdb_path, distance=5.0):
    """
    Returns residue sequence positions within `distance` Å of any ligand in the structure.
    
    Args:
        pdb_path (str): Path to the PDB file.
        distance (float): Distance threshold in Å.

    Returns:
        List[int]: Unique residue numbers of residues near the ligand.
    """
    # Clear current session
    cmd.reinitialize()

    # Load structure
    cmd.load(pdb_path, "protein")

    # Try to select the ligand(s)
    cmd.select("ligand", "organic and not polymer")  # common way to define ligands

    # Select residues within distance
    cmd.select("near_ligand", f"(polymer and (br. all within {distance} of ligand))")

    # Get unique residue identifiers (resi) from selection
    resis = set()
    stored = []
    cmd.iterate("near_ligand", "stored.append(resi)", space={'stored': stored})
    return sorted(set(int(res.strip()) for res in stored if res.strip().isdigit()))

def detect_second_shell(pdb_path, input_residues, distance=5.0):
    """
    Returns all residue sequence numbers (as ints) in `chain` that lie within
    `distance` Å of any residue in `input_residues`, excluding the input set itself.

    Args:
        pdb_path (str): path to .pdb
        input_residues (List[int]): list of residue numbers (e.g. [2,5,7])
        distance (float): cutoff in Å
        chain (str): chain identifier to restrict the search

    Returns:
        List[int]: sorted unique residue numbers in the second shell
    """
    # start fresh
    cmd.reinitialize()
    cmd.load(pdb_path, "protein")

    # turn input into set of strings for selection + keep int set for exclusion
    input_set_str = {str(r) for r in input_residues}
    input_set_int = set(input_residues)

    neighbors_set = set()

    for resi in input_set_str:
        # select that residue in the specified chain
        cmd.select("sel_target", f"polymer and resi {resi}")

        # find all polymer atoms within distance of that selection, still in chain
        cmd.select(
            "sel_neighbors",
            f"polymer and (br. all within {distance} of sel_target)"
        )

        # collect their residue numbers
        stored = []
        cmd.iterate("sel_neighbors", "stored.append(resi)", space={"stored": stored})

        # clean & cast to int (filter out non‑numeric resi like insertion codes)
        for r in stored:
            r_clean = r.strip()
            if r_clean.isdigit():
                neighbors_set.add(int(r_clean))

    # remove the input residues themselves
    result = sorted(neighbors_set - input_set_int)
    return result

def list_to_contig_map(chain_id, seq_length, active_site, missing_residues, start, 
                       segment_extension,n_termini_extension, c_termini_extension, 
                       conservative_RF =False, DSSP_string = None, user_defined_contig_map = None):
    # If you want to have a conservative and an aggresive rf diffusion here's where you do it
    # I might calculate a DSSP string and remove every residue that is in a helix or sheet from the redesign list
    # Account for missing residues in the PDB file
    # Check if the user gave a contig map already
    if user_defined_contig_map is not None:
        print("CONTIG WARNING: Using user defined contig map for RF diffusion redesign.")
        return user_defined_contig_map
    
    # If there are any missing residue positions, we add them to the active site for redesign
    if missing_residues:
        for pos in missing_residues:
            if pos not in active_site and pos > start:
                active_site.append(pos-1)
    
    # Define segments (accounting for conservative/aggresive redesign)
    if conservative_RF:
        segments = define_conservative_segments(active_site, DSSP_string) # TO BE IMPLEMENTED
    else:
        segments = define_aggresive_segments(active_site) 
        print(f"Found segments: {segments}")
    
    # Find chain breaks (TO BE IMPLEMENTED)

    # Assemble contig map
    elements = [rf"[\'{n_termini_extension}-{n_termini_extension}"]  # N-termini extension
    element_start = start
    for seg_start, seg_length in segments:
        element = f"{chain_id}{element_start}-{seg_start}"
        extension = f"{seg_length}-{seg_length + segment_extension}"
        elements.append(element)
        elements.append(extension)
        element_start = seg_start + seg_length + 1
    element = f"{chain_id}{element_start}-{seq_length}"
    elements.append(element)
    if c_termini_extension is None:
        elements[-1] = elements[-1] + rf"\']"
    else:
        elements.append(rf"{c_termini_extension}-{c_termini_extension}\']")
    contig_map = ",".join(elements)
    return contig_map,segments


def define_conservative_segments(active_site, DSSP_string):
    # Define segments based on active site and DSSP string
    return

def define_aggresive_segments(active_site):
    """
    Groups consecutive numbers in active_site into (start, length) tuples.
    """
    # 1. Sort and remove duplicates to ensure logic works
    active_site = sorted(list(set(active_site)))
    
    segments = []
    i = 0
    n = len(active_site)
    
    while i < n:
        segment_start = active_site[i]
        segment_length = 1
        
        # Look ahead: while the next number exists AND is exactly current + 1
        while (i + 1 < n) and (active_site[i + 1] == active_site[i] + 1):
            segment_length += 1
            i += 1 # Advance the index as we consume consecutive numbers
            
        segments.append((segment_start, segment_length))
        
        # Move to the next number to start a new segment check
        i += 1

    return segments




### RF diffusion ###############################################################################################################################    

def run_rf1(
        output_path,
        input_pdb,
        contig_map,
        num_designs,
        T,
        path_to_RF1_script= "~/RFdiffusion/scripts/run_inference.py",
        path_to_RF1_env="/home/eduardo/mambaforge/envs/"
):
    
    
    # Create output directory if it doesn't exist
    Path(f"{output_path}/RF1_designs").mkdir(parents=True, exist_ok=True)

    # Adapt contig map format for RF1 (remove brackets and quotes)
    contig_map = str(contig_map).replace(r"\'","").replace(",","/")

    # Prepare RF1 command
    rf1_command = (
        f'export MKL_THREADING_LAYER=GNU && conda run -p {path_to_RF1_env} python '
        f' {path_to_RF1_script} '
        f'inference.output_prefix="{output_path}/RF1_designs/{str(input_pdb).split("/")[-1].split(".pdb")[0]}_RFdesign" '
        f'inference.input_pdb="{input_pdb}" '
        f'contigmap.contigs="{contig_map}" '
        f'inference.num_designs={num_designs} '
        f'diffuser.T={T}'
    )
    subprocess.run(rf1_command, shell=True)
    return

def run_rfAA(
        output_path,
        input_pdb,
        contig_map,
        num_designs,
        T,
        path_to_RFAA_apptainer ="/home/eduardo/allostery/rf_diffusion_all_atom/rf_se3_diffusion.sif",
        path_to_RFAA_script ="/home/eduardo/allostery/rf_diffusion_all_atom/run_inference.py" ,
        path_to_RFAA_weights = "/home/eduardo/allostery/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt",
        inference_ligand='UNL',
        design_startnum=0,
        deterministic=True
):



    # Create output directory if it doesn't exist
    Path(f"{output_path}/RFallatom_designs").mkdir(parents=True, exist_ok=True)

    # Prepare RF-AA command
    rfAA_command = (
        f'apptainer exec --nv {path_to_RFAA_apptainer} '
        f'python -u {path_to_RFAA_script} '  
        f'inference.output_prefix={output_path}/RFallatom_designs/{str(input_pdb).split("/")[-1].split(".pdb")[0]}_RFdesign '
        f'inference.input_pdb={input_pdb} '
        f'inference.num_designs={num_designs} '
        f'diffuser.T={T} '
        f'contigmap.contigs={contig_map} '
        f'inference.ligand={inference_ligand} '
        f'inference.design_startnum={design_startnum} '
        f'inference.ckpt_path={path_to_RFAA_weights} '
        f'inference.deterministic={str(deterministic).lower()} ' # Hydra/RFAA usually expects lowercase true/false
    )
    print(rfAA_command)
    subprocess.run(rfAA_command, shell=True)
    return

def run_rf3():

    # Skip if output already exists and is not empty
    if os.path.exists(f"{output_path}/RF3_designs") and os.listdir(f"{output_path}/RF3_designs"):
        print(f"RF3 designs already exist in {output_path}/RF3_designs. Skipping RF3 execution.")
        return
    
    return

### LIGAND MPNN #######################################################################################################################################

def generate_redesign_string(pdb_file, original_seq, segments, start_position, chain_id='A'):
    """
    Identifies all positions in an RFdiffusion PDB that are NOT part of the 
    fixed segments (motifs). These positions are targets for LigandMPNN.
    """
    pdb_seq = extract_sequence(str(pdb_file))["coordinate_sequence"]
    
    # 1. Extract the sequences of the fixed motifs from the original
    fixed_motif_seqs = []
    start = start_position
    for seg_start, seg_length in segments:
        # Convert 1-based start to 0-based index
        fixed_motif_seqs.append(original_seq[start-1:seg_start-1])
        start = seg_start + seg_length
    fixed_motif_seqs.append(original_seq[start-1:])


    # 2. Create a mask for the PDB sequence (False = Redesign, True = Fixed)
    # We initialize everything as False (to be redesigned)
    fixed_mask = [False] * len(pdb_seq)
    
    current_search_index = 0
    
    # 3. Find each fixed motif in the PDB sequence sequentially
    print(fixed_motif_seqs)
    for motif in fixed_motif_seqs:
        print(motif)
        # Find the motif starting from the end of the previous found motif
        found_index = pdb_seq.find(motif, current_search_index)
        
        if found_index == -1:
            # Critical Error: RFdiffusion failed to preserve a fixed segment
            raise ValueError(f"Fixed segment '{motif}' lost in PDB sequence!")
        
        # Mark these residues as Fixed (True)
        for i in range(found_index, found_index + len(motif)):
            fixed_mask[i] = True
            
        # Update search start to ensure we look *after* this motif for the next one
        current_search_index = found_index + len(motif)

    # 4. Generate the redesign string for everything NOT masked as fixed
    redesign_list = []
    for i, is_fixed in enumerate(fixed_mask):
        if not is_fixed:
            # Calculate PDB residue number
            res_num = start_position + i
            redesign_list.append(f"{chain_id}{res_num}")
            redesign_string = " ".join(redesign_list)

    return redesign_string

def generate_MPNN_jsons(pdb_folder, output_json,
                                    user_redesign_list, first_shell, chain_id, original_seq,segments, start_position,
                                    first_shell_only=True):
    
    # If output already exists and is not empty skip this function
    if os.path.exists(output_json) and any(os.scandir(output_json)) > 0 :
        print(f"MPNN JSON files already exists at {output_json}. Skipping generation.")
        return
    
    # Create output directory if it doesn't exist
    Path(output_json).mkdir(parents=True, exist_ok=True)
    
    # If first shell, generate the list for first shell
    if first_shell_only:
        first_shell_list = [f'{chain_id}{res}'for res in first_shell]
        first_shell_string = " ".join(first_shell_list)
    
    # Parse the pdb folder and generate json with pdb paths
    pdb_files = list(Path(pdb_folder).glob("*.pdb"))

    data_dict = {}
    redesigned_residues_dict = {}

    for pdb_file in pdb_files:

        # multi pdb
        data_dict[str(pdb_file)] = ""

        # multiredesigns
        if first_shell_only:
            redesigned_residues_dict[str(pdb_file)] = first_shell_list
        else:
            redesigned_residues_dict[str(pdb_file)] = generate_redesign_string(pdb_file, original_seq, segments, start_position, chain_id='A')

    # Generate json files
    json.dump(data_dict, open(f"{output_json}/pdb_paths_multi.json", "w"))
    json.dump(redesigned_residues_dict, open(f"{output_json}/redesigned_residues_multi.json", "w"))
    return

def run_ligand_MPNN(output_folder, num_designs, n_batches, path_to_ligand_MPNN_script, path_to_ligand_MPNN_env,json_path,
                    model_type, temperature, bias_aa_global, omit_aa_global, side_chain_context,model_path):
    
    # Check if output already exists and is not empty
    if os.path.exists(output_folder) and any(f.name.endswith('seqs') for f in os.scandir(output_folder)):
        print(f"MPNN designs already exists at {output_folder}. Skipping generation.")
        return
    
    # Prepare ligand MPNN command
    ligand_MPNN_command = (
        f'conda run -p {path_to_ligand_MPNN_env} python '
        f' {path_to_ligand_MPNN_script} '
        f'--pdb_path_multi {json_path}/pdb_paths_multi.json '
        f'--out_folder {output_folder} '
        f'--save_stats 1 '
        f'--batch_size {num_designs} '
        f'--number_of_batches {n_batches} '
        f'--model_type {model_type} '
        f'--checkpoint_ligand_mpnn {model_path} '
        f'--temperature {temperature} '
        f'--ligand_mpnn_use_side_chain_context {side_chain_context} '
        f'--redesigned_residues_multi "{json_path}/redesigned_residues_multi.json" '
    )

    # Append global bias if given
    if omit_aa_global:
        ligand_MPNN_command += f'--omit_AA {omit_aa_global} '

    # Append global bias if given
    if bias_aa_global:
        ligand_MPNN_command += f'--bias_AA {bias_aa_global}  '

    # Append fixed residues if the file exists
    fixed_residues_file = os.path.join(json_path, "fix_residues_multi.json")
    if os.path.exists(fixed_residues_file):
        ligand_MPNN_command += f'--fixed_residues_multi "{fixed_residues_file}" '
    
    # Append multi omit AA per residue if the file exists
    omit_AA_per_residue_multi_file = os.path.join(json_path, "omit_AA_per_residue_multi.json")
    if os.path.exists(omit_AA_per_residue_multi_file):
        ligand_MPNN_command += f'--omit_AA_per_residue_multi "{omit_AA_per_residue_multi_file}" '
    
    # Append multi bias AA per residue if the file exists
    bias_AA_per_residue_multi_file = os.path.join(json_path, "bias_AA_per_residue_multi.json")
    if os.path.exists(bias_AA_per_residue_multi_file):
        ligand_MPNN_command += f'--bias_AA_per_residue_multi "{bias_AA_per_residue_multi_file}" '

    subprocess.run(ligand_MPNN_command, shell=True)

    return

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

def run_ESMfold(input_csv,output_folder,path_to_ESM_env,path_to_ESM_script):
    # 1. Convert to Path object first
    out_path = Path(output_folder)

    # 2. Check using .is_dir() (safer than .exists) and .iterdir()
    if out_path.is_dir() and any(out_path.iterdir()):
        print(f"ESM predictions already exists at {output_folder}. Skipping generation.")
        return
    
    # Prepare command for running the ESM script
    ESM_command = (
        f'conda run -p {path_to_ESM_env} python3 -u {path_to_ESM_script} --input_csv {input_csv} --output_folder {output_folder}'
    )
    subprocess.run(ESM_command, shell=True, check=True)
    return

def compute_msa():
    return

def AF_json_generator_wo_MSA(row, json_path):
    file_name = row['file_ID']
    seq = row['sequence']

    json_data = {
        "name": f"{file_name}",
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq,
                    "unpairedMsa": f">dummy\n{seq}\n",  
                    "pairedMsa": "",                    
                    "templates": []
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }

    with open(f"{json_path}/temp.json", "w") as f:
        json.dump(json_data, f, indent=2)

def AF_json_generator_w_MSA(row, json_path):
    return

def AF_json_generator_cofold(row,ligand_SMILES, json_path):
    file_name = row['file_ID']
    seq = row['sequence']
    partial_T = row['partial_T']

    json_data = {
        "name": f"{file_name}",
        "sequences": [
            {
                "protein": {
                    "id": "A",
                    "sequence": seq,
                    "unpairedMsa": f">dummy\n{seq}\n",  
                    "pairedMsa": "",                    
                    "templates": []
                }
            },
            {
                "ligand": {
                    "id": "L",
                    "smiles": ligand_SMILES
                }
            }
        ],
        "modelSeeds": [1],
        "dialect": "alphafold3",
        "version": 1
    }

    with open(f"{json_path}/temp.json", "w") as f:
        json.dump(json_data, f, indent=2)
    return

def run_AlphaFold3(row, output_path,msa_flag=False):
    # Create temporary JSON file
    if msa_flag:
        AF_json_generator_w_MSA(row, output_path)
    else:
        AF_json_generator_wo_MSA(row, output_path)
    # Run AF3
    af_command = (
    f"conda run -p /home/ingemar/anaconda3/envs/alphafold3 python /mnt/data/alphafold3/run_alphafold.py "
    f"--json_path={output_path}/temp.json "
    f"--output_dir={output_path} "
    f"--db_dir=/mnt/data/alphafold3/alphafold_databases/ "
    f"--model_dir=/mnt/data/alphafold3/models/ "
    f"--num_diffusion_samples=1"
    )
    af_process = subprocess.run(af_command, shell = True, text=True)

def run_chai(row, output_path):
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
    chai_command =f" conda run -p /home/ingemar/miniforge3/envs/chai chai-lab fold {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return 

def run_chai_w_MSA(row, output_path):
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
    chai_command =f" conda run -p /home/ingemar/miniforge3/envs/chai chai-lab fold --use-msa-server --use-templates-server {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return 

def run_chai_cofold(row, output_path,ligand_smiles):
    #
    Path(f"{output_path}_cofold").mkdir(parents=True, exist_ok=True)
    #
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}_cofold/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
        f.write(f">ligand|{file_name}\n")      # Header line
        f.write(f"{ligand_smiles}\n")         # The actual ligand SMILES
    # Run chai
    chai_command =f" conda run -p /home/ingemar/miniforge3/envs/chai chai-lab fold {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return 

def run_chai_cofold_w_MSA(row, output_path,ligand_smiles):
    #
    Path(f"{output_path}_cofold").mkdir(parents=True, exist_ok=True)
    #
    sequence = row['sequence']
    file_name = row['file_ID']
    # Create the full path for the temporary FASTA file
    fasta_path = f"{output_path}_cofold/temp.fa"

    # Write the sequence to the FASTA file
    with open(fasta_path, 'w') as f:
        f.write(f">protein|{file_name}\n")      # Header line
        f.write(f"{sequence}\n")         # The actual sequence
        f.write(f">ligand|{file_name}\n")      # Header line
        f.write(f"{ligand_smiles}\n")         # The actual ligand SMILES
    # Run chai
    chai_command =f" conda run -p /home/ingemar/miniforge3/envs/chai chai-lab fold --use-msa-server --use-templates-server {fasta_path} {output_path}/{file_name.lower()}" 
    chai_process = subprocess.run(chai_command, shell=True)

    return 

def second_prediction_round(df,model_flag,MSA_flag,ligand_smiles, output_path):
    if model_flag not in ['AF3','CHAI']:
        print("Error: Model flag must be either 'AF3' or 'CHAI'")
        sys.exit(1)
    for index, row in df.iterrows():
        if model_flag == 'AF3':
            if MSA_flag:
                run_AlphaFold3(row, output_path, msa_flag=True)
            else:
                run_AlphaFold3(row, output_path, msa_flag=False)
        elif model_flag == 'CHAI':
            if ligand_smiles:
                if MSA_flag:
                    run_chai_cofold_w_MSA(row, output_path, ligand_smiles)
                else:
                    run_chai_cofold(row, output_path, ligand_smiles)
            else:
                if MSA_flag:
                    run_chai_w_MSA(row, output_path)
                else:
                    run_chai(row, output_path)
                    
    if model_flag=="AF3":
        process_AF3_folder(output_path)
    if model_flag=="CHAI":
        process_chai_folder(output_path)
    return
### PROTEIN METRICS #########################################################################################################################################


def extract_plddt(file_path):
    """
    Extracts the mean and standard deviation of pLDDT scores 
    from a structure file using Biotite.
    """
    try:
        # FIX: Explicitly ask for 'b_factor' using extra_fields
        struct = bsio.load_structure(file_path, extra_fields=["b_factor"])
        
        # Check if b_factor data was successfully loaded
        if not hasattr(struct, "b_factor"):
            # Fallback: Sometimes PDB parsers attach it to 'occupancy' by mistake 
            # if columns are misaligned, but usually it's just missing from the load.
            print(f"Warning: 'b_factor' column could not be loaded from {file_path}")
            return 0.0, 0.0
            
        # Extract the array of pLDDT scores
        plddt_values = struct.b_factor
        
        # Return Mean and Stdev
        return np.mean(plddt_values), np.std(plddt_values)

    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return 0.0, 0.0



# Define the dictionary at the module level for efficiency
THREE_TO_ONE_MAP = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def three_to_one(res_names):
    """
    Converts a list or array of 3-letter residue codes to a single 1-letter string.
    Handles non-standard residues by mapping them to 'X'.
    """
    # Use .get() to return 'X' for any residue not in the dictionary
    return "".join([THREE_TO_ONE_MAP.get(res, 'X') for res in res_names])

def get_coords(pdb_path):
    """
    Helper to get (N, 3) coordinates of CA atoms and the sequence string.
    Ensures coordinates are float64 and sequence is a string.
    """
    struct = load_structure(pdb_path)
    
    # Filter for Alpha Carbons only
    ca = struct[struct.atom_name == "CA"]
    
    # Explicitly cast to float64 for tm_align compatibility
    coords_64 = ca.coord.astype(np.float64)
    
    # Convert 3-letter code array to single string
    seq_str = three_to_one(ca.res_name)
    
    return coords_64, seq_str

def batch_tm_align(target_path, design_path):
    # 1. Load Target Coords
    target_coords, target_seq = get_coords(target_path)
    
    results = []
    
    mobile_coords, mobile_seq = get_coords(design_path)
    
    # 2. Run TM-align
    # output is a specialized object containing score and rotation matrix
    res = tm_align(mobile_coords, target_coords, mobile_seq, target_seq)
        
    return res.rmsd, res.tm_norm_chain1
def detect_clashes(file_path, clash_distance=2.0, bond_distance=1.2):
    return 0, 0.0
def detect_clashes_2(file_path, clash_distance=2.0, bond_distance=1.2):
    """
    Detects steric clashes using Biotite (Vectorized/Fast).
    Corrected to handle NumPy adjacency matrices.
    """
    try:
        # 1. Load Structure
        structure = load_structure(file_path)
        
        # 2. Filter out Hydrogens (Vectorized)
        mask = (structure.element != "H")
        atoms = structure[mask]
        
        if len(atoms) == 0:
            return 0, 0.0

        # 3. Efficient Neighbor Search
        cell_list = struc.CellList(atoms, cell_size=clash_distance)
        
        # Returns a boolean NumPy array (N x N)
        adjacency = cell_list.create_adjacency_matrix(clash_distance)
        
        # FIX: Use np.where to get indices from a dense NumPy array
        # rows and cols will be arrays of indices where an interaction exists
        rows, cols = np.where(adjacency)
        
        # 4. Filter Pairs
        # Filter for unique pairs (row < col) to remove duplicates and self-loops
        pair_mask = (rows < cols)
        rows = rows[pair_mask]
        cols = cols[pair_mask]
        
        # 5. Calculate Exact Distances & Apply Bond Filter
        coords_A = atoms.coord[rows]
        coords_B = atoms.coord[cols]
        
        # Vectorized Euclidean distance calculation
        diffs = coords_A - coords_B
        dists_sq = np.sum(diffs**2, axis=1)
        dists = np.sqrt(dists_sq)
        
        # Logic: Clash if dist < clash_distance AND dist > bond_distance
        real_clashes_mask = (dists > bond_distance)
        
        # 6. Count Results
        num_clashes = np.sum(real_clashes_mask)
        num_atoms = len(atoms)
        clashes_per_atom = num_clashes / num_atoms if num_atoms > 0 else 0.0
        
        return num_clashes, clashes_per_atom

    except Exception as e:
        print(f"Error in clash detection for {file_path}: {e}")
        return 0, 0.0

def compute_SAPscore():
    return

### BINDING AND POCKET METRICS ###############################################################################################################################

def gnina_minimize_autobox(pdb_file,
                           ligand,
                           output_folder,
                           gnina_path,
                           cnn,
                           exhaustiveness, 
                           autobox_add):
    """
    This function both creates poses for the ligand for ESM predictions as well as computes the gnina scores
    for each structure in the pdb_folder. Returns a dataframe with GNina metrics
    
    :param pdb_folder: Description
    :param ligand: Description
    :param output_folder: Description
    :param gnina_path: Description
    :param cnn: Description
    :param exhaustiveness: Description
    :param autobox_add: Description
    """
    # Create output directory if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)

    # Prepare gnina command
    # Note: Added check to ensure pdb_file is a Path object or string
    pdb_path = Path(pdb_file)
    output_path = Path(output_folder) / f"{pdb_path.stem}_ligand.pdb"

    gnina_command = (
        f'{gnina_path} -r {pdb_file} -l {ligand} --autobox_ligand {pdb_file} '
        f'-o {output_path} --minimize '
        f'--exhaustiveness {exhaustiveness} --cnn {cnn} --autobox_add {autobox_add}'
    )
    
    # Run the command
    gnina_result = subprocess.run(gnina_command, shell=True, capture_output=True, text=True)
    
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
### 3D FILTERING #######################################################################################################################################

def sliding_window_1D_filter(seq, window_size, threshold, aa='A'):
    if len(seq) < window_size:
        return seq.count(aa) >= threshold
    for i in range(len(seq) - window_size + 1):
        window = seq[i : i + window_size]
        if window.count(aa) >= threshold:
            return True 
    return False

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

def threed_params_1_df(folder, output_folder, original_path, clash_distance=2.0, bond_distance=1.2,
                       ligand_path=None, gnina_path="gnina", cnn="default",
                       exhaustiveness=8, autobox_add=4):
    """
    Given a folder with structures, returns the following df: "file_ID", "sequence", "pLDDT mean", "pLDDT std dev", "RMSD", "TMscore", "clashes"
    """
    # Skip if output already exists and is not empty
    if os.path.exists(f"{output_folder}/ESM_filtered_df.csv") and os.path.getsize(f"{output_folder}/ESM_filtered_df.csv") > 0:
        print(f"3D metrics already exists at {output_folder}/ESM_filtered_df.csv. Skipping generation.")
        df = pd.read_csv(f"{output_folder}/ESM_filtered_df.csv")
        return df
    # Initialise df
    elements = ["file_ID", "sequence", "pLDDT_mean", "pLDDT_std", "RMSD", "TMscore", "num_clashes","clashes_per_atom","Gnina_CNNscore","Gnina_affinity"]
    threed_df = pd.DataFrame(columns = elements)

    # Parse folder
    folder = Path(folder)
    for file in folder.iterdir():
        file_path = str(file)
        if not file_path.endswith(".json"):
            # pLDDT
            pLDDT_mean,pLDDT_stdev= extract_plddt(file_path)
            # RMSD
            RMSD, TMscore = batch_tm_align(original_path,file_path)
            # Clashes
            num_clashes,clash_atom = detect_clashes(file_path, clash_distance=clash_distance, bond_distance=bond_distance)
            # Sequence
            seq = extract_sequence(file_path, chain_id='A')["coordinate_sequence"]
            # GNINA scores
            Gnina_CNNscore, Gnina_affinity,vina_affinity,vina_affinity_2 = gnina_minimize_autobox(pdb_file=file,
                                                                    ligand=ligand_path,
                                                                    output_folder=f"{output_folder}/ESM_Gnina_minimized/",
                                                                    gnina_path=gnina_path,
                                                                    cnn=cnn,
                                                                    exhaustiveness=exhaustiveness,
                                                                    autobox_add=autobox_add)
            # Generate row
            row = pd.Series([file.name, seq, pLDDT_mean, pLDDT_stdev, RMSD, TMscore, num_clashes,clash_atom,Gnina_CNNscore,Gnina_affinity], index=list(threed_df.columns))
            threed_df.loc[len(threed_df)] = row

    print(threed_df.head(5))
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


### AUXILIARY FUNCTIONS ###############################################################################################################################
def move_pdbs_to_folder(input_folder, output_folder):
    """
    Moves all PDB files from pdb_folder to output_folder.
    """
    # Create output folder if it doesn't exist
    Path(output_folder).mkdir(parents=True, exist_ok=True)
    
    # Move each PDB file
    for pdb_file in Path(input_folder).glob("*.pdb"):
        new_path = Path(output_folder) / pdb_file.name
        pdb_file.rename(new_path)

def process_chai_folder(folder):
    return

def process_AF3_folder(folder):
    return