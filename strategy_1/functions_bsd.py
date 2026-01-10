"""
Functions for running the script binding_site_designer.py
"""

import warnings
import os
import json
from Bio import PDB
from Bio.PDB import PDBParser, DSSP
from Bio.SeqUtils import seq1
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
        "pdb_id": pdb_path.split('/')[-1].split('.')[0],
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
    return contig_map


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

def list_to_MPNN_json():

    return

def json_generator(data_dict, output_dir,phase_change):
    """
    Write one .jsonl file per entry in the dictionary.

    Each file will be named <key>.jsonl and contain a single line like:
    {"<key>": {"A": [list of values]}}
    """
    os.makedirs(output_dir, exist_ok=True)

    for key, value in data_dict.items():
        json_data = {key.split(".")[0]: {"A": [x - phase_change for x in value]}}
        output_path = os.path.join(output_dir, f"{key}.jsonl")
        with open(output_path, "w") as f:
            f.write(json.dumps(json_data) + "\n")

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

### Folding models ###############################################################################################################################

### 3D FILTERING #######################################################################################################################################

### BINDING AND POCKET METRICS ###############################################################################################################################

### AUXILIARY FUNCTIONS ###############################################################################################################################
