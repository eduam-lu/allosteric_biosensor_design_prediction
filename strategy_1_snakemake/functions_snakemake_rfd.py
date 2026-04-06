"""
Functions for RF Diffusion Binding Site Designer - Snakemake Pipeline

This module provides functions for:
- PDB structure parsing and information extraction
- Active site detection and definition
- Contig map generation for RF Diffusion
- RF Diffusion model execution (RF1, RF Diffusion All-Atom, RF3)
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
import shutil
import pandas as pd
import re
from pymol import cmd 
from pathlib import Path
import subprocess
import math
from scipy import stats
from snakemake.shell import shell


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def load_structure(pdb_path):
    """
    Load PDB or CIF structure files using the appropriate parser.
    
    Args:
        pdb_path (str): Path to PDB or CIF file
        
    Returns:
        Bio.PDB.Structure: Parsed structure object
    """
    pdb_path_str = str(pdb_path).lower()
    
    if pdb_path_str.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    
    structure = parser.get_structure("structure", pdb_path)
    return structure


# ============================================================================
# PDB INFORMATION EXTRACTION
# ============================================================================

def extract_pdb_info(pdb_path, 
                    first_shell_distance=5.0, second_shell_distance=5.0, 
                    user_defined_active_site=None, user_defined_residues=None):
    """
    Extract comprehensive PDB information including sequence, secondary structure, and active site.
    
    Args:
        pdb_path (str): Path to PDB file
        first_shell_distance (float): Distance threshold for first shell detection in Ångström (default 5.0)
        second_shell_distance (float): Distance threshold for second shell detection in Ångström (default 5.0)
        user_defined_active_site (list): User-defined active site residues (optional)
        user_defined_residues (list): User-defined residues to include (optional)
        
    Returns:
        tuple: (pdb_info_dict, printable_pdb_info_dict)
               Both dictionaries contain PDB ID, chain ID, sequences, secondary structure, active site info
    """
    # Extract chain id
    chain_id = get_chain_id(pdb_path)
    # Extract sequence
    sequence_dict = extract_sequence(pdb_path, chain_id)
    # Extract DSSP
    dssp_string = extract_dssp_string(pdb_path)
    # Define active site
    active_site_dict = define_active_site(pdb_path, 
                                          first_shell_distance=first_shell_distance, 
                                          second_shell_distance=second_shell_distance, 
                                          user_defined_active_site=user_defined_active_site, 
                                          user_defined_residues=user_defined_residues)
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
    Get the chain ID from a PDB file.
    
    If multiple chains exist, returns the first one found.
    
    Args:
        pdb_path (str): Path to PDB file
        
    Returns:
        str: Chain ID (e.g., 'A', 'B'), or None if not found
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return None

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
    Extract sequence information from a PDB or CIF file.
    
    Combines SEQRES (if available) and coordinate-based sequences,
    handling missing residues and gaps appropriately.
    
    Args:
        pdb_path (str): Path to PDB or CIF file
        chain_id (str): Chain ID to extract (default 'A')
        
    Returns:
        dict: Dictionary containing:
            - pdb_id: PDB identifier
            - seqres: SEQRES sequence from header
            - coordinate_sequence: Sequence from coordinates
            - residue_count: Number of residues
            - found_residues: Number of non-missing residues
            - missing_residues: Number of missing residues
            - missing_positions: List of missing positions
            - start_residue_number: Position of first residue
            - sequence_length: Total sequence length
    """
    # 1. Select the correct parser
    if str(pdb_path).endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)

    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {"error": f"Structure Parsing Error: {str(e)}"}

    # 2. Get SEQRES sequence
    try:
        seqres = get_seqres_sequence(structure, chain_id)
    except Exception:
        seqres = None

    # 3. Get Coordinate Sequence
    try:
        if chain_id in structure[0]:
            coord_seq, res_count = get_coordinate_sequence(structure[0][chain_id])
        else:
            coord_seq, res_count = ("", 0)
    except KeyError:
        return {"error": f"Chain {chain_id} not found in structure"}

    # 4. Fallback Logic
    if not seqres:
        seqres = coord_seq  # Fallback to coordinate sequence if SEQRES not found

    # 5. Alignment / Padding Logic
    if len(coord_seq) < len(seqres):
        coord_seq += '-' * (len(seqres) - len(coord_seq))
    
    missing_residues = coord_seq.count('-')
    found_residues = len(coord_seq) - missing_residues
    missing_positions = [i+1 for i, res in enumerate(coord_seq) if res == '-']

    # 6. Find start position
    first_residue_pos = next((i for i, res in enumerate(coord_seq) if res != '-'), None)
    if first_residue_pos is not None:
        start_residue_number = first_residue_pos + 1  # 1-based indexing
    else:
        start_residue_number = None 

    return {
        "pdb_id": os.path.basename(pdb_path).split('.')[0],
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
    Retrieve the SEQRES sequence from the PDB header.
    
    Args:
        structure (Bio.PDB.Structure): Parsed structure object
        chain_id (str): Chain ID to extract (default 'A')
        
    Returns:
        str: SEQRES sequence string, or None if not found
    """
    if 'seqres' in structure.header:
        seqres_dict = structure.header['seqres']
        if chain_id in seqres_dict:
            return str(seqres_dict[chain_id])
    return None


def get_coordinate_sequence(chain):
    """
    Extract amino acid sequence from chain coordinates.
    
    Handles missing residues by inserting gaps ('-') and accounts for
    N-terminal gaps by padding from position 1.
    
    Args:
        chain (Bio.PDB.Chain): Chain object
        
    Returns:
        tuple: (sequence_string, residue_count)
               sequence_string: 1-letter codes with '-' for gaps
               residue_count: Number of actual residues found
    """
    # Filter for standard amino acids (ignore waters/ligands)
    residues = [res for res in chain if PDB.is_aa(res)]
    
    if not residues:
        return "", 0

    # Extract residue numbers (id[1]) and 1-letter codes
    res_map = {res.id[1]: seq1(res.resname) for res in residues}
    
    start_res_num = min(res_map.keys())
    end_res_num = max(res_map.keys())
    
    # Iterate from 1 (standard start) to the last observed residue
    # This automatically handles N-term gaps (1 to start) and internal gaps
    sequence_parts = []
    for i in range(1, end_res_num + 1):
        if i in res_map:
            sequence_parts.append(res_map[i])
        else:
            sequence_parts.append('-')  # Gap detected
            
    final_seq = "".join(sequence_parts)
    residue_count = len(residues)  # Count of actual physical residues
    
    return final_seq, residue_count


def extract_dssp_string(input_file):
    """
    Extract secondary structure assignments using DSSP.
    
    Args:
        input_file (str): Path to PDB file
        
    Returns:
        str: DSSP secondary structure string (H=helix, E=sheet, C=coil, etc.)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_file)
    model = structure[0]  # DSSP works on one model

    dssp = DSSP(model, input_file)
    dssp_codes = ''.join([dssp[key][2] for key in dssp.keys()])

    return dssp_codes


def define_active_site(pdb_path, user_defined_active_site=None, user_defined_residues=None, 
                      first_shell_distance=5.0, second_shell_distance=5.0):
    """
    Define the active site based on proximity to ligand.
    
    Detects residues close to any ligand in the structure, optionally
    extending to a second shell and including user-defined residues.
    
    Args:
        pdb_path (str): Path to PDB file
        user_defined_active_site (list): Pre-defined active site residues (optional)
        user_defined_residues (list): Additional residues to include (optional)
        first_shell_distance (float): Distance for first shell in Ångström (default 5.0)
        second_shell_distance (float): Distance for second shell in Ångström (default 5.0)
        
    Returns:
        dict: Dictionary with keys:
            - active_site: List of active site residue numbers
            - first_shell: First shell residues
            - second_shell: Second shell residues
            - user_defined_active_site: User-provided active site
            - user_defined_residues: User-provided additional residues
    """
    # If active site already defined, use it and add user residues
    if user_defined_active_site is not None:
        active_site = user_defined_active_site
        if user_defined_residues:
            for i in user_defined_residues:
                if i not in active_site:
                    active_site.append(i)
    else:
        # Auto-detect based on ligand proximity
        first_shell = detect_first_shell(pdb_path, distance=first_shell_distance)
        second_shell = detect_second_shell(pdb_path, first_shell, distance=second_shell_distance)
        active_site = first_shell

        # Include user-defined residues
        if user_defined_residues is not None:
            for i in user_defined_residues:
                if i not in active_site:
                    active_site.append(i)

    first_shell = detect_first_shell(pdb_path, distance=first_shell_distance)
    second_shell = detect_second_shell(pdb_path, first_shell, distance=second_shell_distance)

    return {
        "active_site": active_site,
        "first_shell": first_shell,
        "second_shell": second_shell,
        "user_defined_active_site": user_defined_active_site,
        "user_defined_residues": user_defined_residues
    }


def detect_first_shell(pdb_path, distance=5.0):
    """
    Detect residues within specified distance of any ligand.
    
    Uses PyMOL to select ligand atoms and find nearby protein residues.
    
    Args:
        pdb_path (str): Path to PDB file
        distance (float): Distance threshold in Ångström (default 5.0)

    Returns:
        list: Sorted list of residue numbers within the specified distance
    """
    cmd.reinitialize()
    cmd.load(pdb_path, "protein")
    cmd.select("ligand", "organic and not polymer")
    cmd.select("near_ligand", f"(polymer and (br. all within {distance} of ligand))")

    stored = []
    cmd.iterate("near_ligand", "stored.append(resi)", space={'stored': stored})
    return sorted(set(int(res.strip()) for res in stored if res.strip().isdigit()))


def detect_second_shell(pdb_path, input_residues, distance=5.0):
    """
    Detect residues within distance of first shell residues.
    
    Returns all residues near the input residue set, excluding the input set itself.
    
    Args:
        pdb_path (str): Path to PDB file
        input_residues (list): First shell residue numbers
        distance (float): Distance threshold in Ångström (default 5.0)

    Returns:
        list: Sorted list of second shell residue numbers
    """
    cmd.reinitialize()
    cmd.load(pdb_path, "protein")

    input_set_str = {str(r) for r in input_residues}
    input_set_int = set(input_residues)
    neighbors_set = set()

    for resi in input_set_str:
        cmd.select("sel_target", f"polymer and resi {resi}")
        cmd.select("sel_neighbors", f"polymer and (br. all within {distance} of sel_target)")

        stored = []
        cmd.iterate("sel_neighbors", "stored.append(resi)", space={"stored": stored})

        for r in stored:
            r_clean = r.strip()
            if r_clean.isdigit():
                neighbors_set.add(int(r_clean))

    result = sorted(neighbors_set - input_set_int)
    return result


# ============================================================================
# CONTIG MAP GENERATION
# ============================================================================

def list_to_contig_map(chain_id, seq_length, active_site, missing_residues, start, 
                       segment_extension, n_termini_extension, c_termini_extension, 
                       conservative_RF=False, DSSP_string=None, user_defined_contig_map=None):
    """
    Generate contig map string for RF Diffusion from active site definition.
    
    Converts residue list into RF Diffusion contig format specifying
    which regions should be designed vs fixed.
    
    Args:
        chain_id (str): Protein chain ID
        seq_length (int): Total sequence length
        active_site (list): Active site residue numbers
        missing_residues (list): Missing residue positions
        start (int): Starting residue number
        segment_extension (int): Extension around active site segments
        n_termini_extension (int): N-terminal extension length
        c_termini_extension (int): C-terminal extension length (optional)
        conservative_RF (bool): Use conservative redesign mode (default False)
        DSSP_string (str): Secondary structure string (optional)
        user_defined_contig_map (str): Pre-defined contig map (optional)
        
    Returns:
        tuple: (contig_map_string, segments)
               contig_map_string: RF Diffusion format contig specification
               segments: List of (start, length) tuples for designed segments
    """
    # If user-defined contig map provided, use it
    if user_defined_contig_map is not None:
        print("CONTIG: Using user-defined contig map for RF Diffusion redesign.")
        return user_defined_contig_map, None
    
    # If there are missing residues, add them to active site for redesign
    if missing_residues:
        for pos in missing_residues:
            if pos not in active_site and pos > start:
                active_site.append(pos - 1)
    
    # Define segments based on redesign strategy
    if conservative_RF:
        segments = define_conservative_segments(active_site, DSSP_string)
    else:
        segments = define_aggresive_segments(active_site)
        print(f"Designed segments: {segments}")
    
    # Assemble contig map in RF Diffusion format
    elements = [rf"[\'{n_termini_extension}-{n_termini_extension}"]
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
    return contig_map, segments


def define_conservative_segments(active_site, DSSP_string):
    """
    Define redesign segments using conservative strategy with secondary structure.
    
    Filters active site residues to exclude secondary structure elements.
    
    Args:
        active_site (list): Active site residue numbers
        DSSP_string (str): Secondary structure assignment string
        
    Returns:
        list: Filtered list of (start, length) segment tuples
    """
    # Conservative strategy: exclude residues in defined secondary structures
    return define_aggresive_segments(active_site)


def define_aggresive_segments(active_site):
    """
    Group consecutive residues into segments.
    
    Converts an active site residue list into (start, length) tuples
    representing contiguous segments for redesign.
    
    Args:
        active_site (list): Residue numbers to group
        
    Returns:
        list: List of (segment_start, segment_length) tuples
    """
    active_site = sorted(list(set(active_site)))
    
    segments = []
    i = 0
    n = len(active_site)
    
    while i < n:
        segment_start = active_site[i]
        segment_length = 1
        
        # Look ahead for consecutive numbers
        while (i + 1 < n) and (active_site[i + 1] == active_site[i] + 1):
            segment_length += 1
            i += 1
            
        segments.append((segment_start, segment_length))
        i += 1

    return segments


# ============================================================================
# RF DIFFUSION EXECUTION  
# ============================================================================

def run_rf1(output_path, input_pdb, contig_map, num_designs, T,
            path_to_RF1_script="~/RFdiffusion/scripts/run_inference.py",
            path_to_RF1_env="/home/eduardo/mambaforge/envs/"):
    """
    Execute RF Diffusion (original RF1 model).
    
    Args:
        output_path (str): Output directory for designed structures
        input_pdb (str): Path to input PDB file
        contig_map (str): Contigs specifying design regions
        num_designs (int): Number of designs to generate
        T (int): Diffusion temperature/steps  
        path_to_RF1_script (str): Path to RF1 inference script
        path_to_RF1_env (str): Path to RF1 conda environment
    """
    Path(f"{output_path}/RF1_designs").mkdir(parents=True, exist_ok=True)

    contig_map = str(contig_map).replace(r"\'", "").replace(",", "/")

    rf1_command = (
        f'export MKL_THREADING_LAYER=GNU && conda run -p {path_to_RF1_env} python '
        f'{path_to_RF1_script} '
        f'inference.output_prefix="{output_path}/RF1_designs/{str(input_pdb).split("/")[-1].split(".pdb")[0]}_RFdesign" '
        f'inference.input_pdb="{input_pdb}" '
        f'contigmap.contigs="{contig_map}" '
        f'inference.num_designs={num_designs} '
        f'diffuser.T={T}'
    )
    subprocess.run(rf1_command, shell=True)


def run_rfAA(output_path, input_pdb, contig_map, num_designs, T,
             path_to_RFAA_apptainer="/home/eduardo/allostery/rf_diffusion_all_atom/rf_se3_diffusion.sif",
             path_to_RFAA_script="/home/eduardo/allostery/rf_diffusion_all_atom/run_inference.py",
             path_to_RFAA_weights="/home/eduardo/allostery/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt",
             inference_ligand='UNL', design_startnum=0, deterministic=True):
    """
    Execute RF Diffusion All-Atom model.
    
    Supports ligand-aware all-atom protein design with specified ligand residue.
    
    Args:
        output_path (str): Output directory for designed structures
        input_pdb (str): Path to input PDB file with ligand
        contig_map (str): Contigs specifying design regions
        num_designs (int): Number of designs to generate
        T (int): Diffusion temperature/steps
        path_to_RFAA_apptainer (str): Path to RF-AA Apptainer/Singularity image
        path_to_RFAA_script (str): Path to RF-AA inference script
        path_to_RFAA_weights (str): Path to RF-AA model weights
        inference_ligand (str): Ligand residue name for design awareness
        design_startnum (int): Starting design number (default 0)
        deterministic (bool): Use deterministic mode (default True)
    """
    Path(f"{output_path}/RFallatom_designs").mkdir(parents=True, exist_ok=True)

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
        f'inference.deterministic={str(deterministic).lower()}'
    )
    print(rfAA_command)
    subprocess.run(rfAA_command, shell=True)


def generate_rf3_yaml(output_yaml_path, input_pdb, contig_map, ligand_resname, num_designs_batch, 
                      chain_id='A', checkpoint_path='rfd3', batch_size=1, 
                      use_classifier_free_guidance=False, cfg_scale=1.5, num_timesteps=200, 
                      step_scale=1.5, noise_scale=1.003, gamma_0=0.6, gamma_min=1.0, 
                      dump_trajectories=False, prevalidate_inputs=False, low_memory_mode=False):
    """
    Generate RF Diffusion 3 configuration YAML file.
    
    Args:
        output_yaml_path (str): Path to save the YAML configuration
        input_pdb (str): Path to input PDB structure
        contig_map (str): Contig string specifying design regions
        ligand_resname (str): Ligand residue name (e.g., 'UNL', 'LIG')
        num_designs_batch (int): Number of designs per batch
        chain_id (str): Protein chain ID (default 'A')
        checkpoint_path (str): Path to RF3 model checkpoint (default 'rfd3')
        batch_size (int): Diffusion batch size (default 1)
        use_classifier_free_guidance (bool): Enable CFG (default False)
        cfg_scale (float): CFG guidance scale (default 1.5) 
        num_timesteps (int): Number of diffusion timesteps (default 200)
        step_scale (float): Step scaling factor (default 1.5)
        noise_scale (float): Noise scaling factor (default 1.003)
        gamma_0 (float): Initial gamma parameter (default 0.6)
        gamma_min (float): Minimum gamma parameter (default 1.0)
        dump_trajectories (bool): Save diffusion trajectories (default False)
        prevalidate_inputs (bool): Validate inputs before design (default False)
        low_memory_mode (bool): Enable low memory consumption mode (default False)
    """
    def _clean_str(value):
        """Strip outer list brackets, backslashes, and quotes from a string."""
        if isinstance(value, (list, tuple)):
            value = value[0] if len(value) == 1 else ",".join(map(str, value))
        s = str(value).strip()
        if s.startswith("[") and s.endswith("]"):
            s = s[1:-1].strip()
        s = re.sub(r'^[\\\'\" ]+|[\\\'\" ]+$', '', s)
        return s

    normalized_contig = _clean_str(contig_map)
    normalized_ligand = "UNL"

    # Build RF3 inputs YAML
    input_specs = {
        "spec-1": {
            "input": input_pdb,
            "contig": normalized_contig,
            "select_fixed_atoms": {
                "UNL": ""
            },
            "ligand": normalized_ligand
        }
    }

    with open(output_yaml_path, 'w') as f:
        yaml.safe_dump(input_specs, f, sort_keys=False)
    
    print(f"RF3 configuration YAML generated: {output_yaml_path}")


def run_rfd3(output_path, input_pdb, contig_map, pdb_info, num_designs, chain_id, 
             path_to_RF3_env=None, ligand_resname='LIG', checkpoint_path='rfd3', 
             batch_size=1, use_classifier_free_guidance=False, cfg_scale=1.5, 
             num_timesteps=200, step_scale=1.5, noise_scale=1.003, gamma_0=0.6, 
             gamma_min=1.0, dump_trajectories=False, prevalidate_inputs=False, 
             low_memory_mode=False):
    """
    Execute RF Diffusion 3 design pipeline.
    
    Latest RF Diffusion model with enhanced capabilities including classifier-free guidance.
    
    Args:
        output_path (str): Output directory for designed structures
        input_pdb (str): Path to input PDB structure
        contig_map (str): Contig string specifying design regions
        pdb_info (dict): Dictionary with PDB information (required for active site)
        num_designs (int): Number of designs to generate
        chain_id (str): Protein chain ID
        path_to_RF3_env (str): Path to RF3 environment (conda prefix or path)
        ligand_resname (str): Ligand residue name for design awareness (default 'LIG')
        checkpoint_path (str): Path to RF3 model checkpoint (default 'rfd3')
        batch_size (int): Diffusion batch size (default 1)
        use_classifier_free_guidance (bool): Enable CFG (default False)
        cfg_scale (float): CFG guidance scale (default 1.5)
        num_timesteps (int): Number of diffusion timesteps (default 200)
        step_scale (float): Step scaling factor (default 1.5)
        noise_scale (float): Noise scaling factor (default 1.003)
        gamma_0 (float): Initial gamma parameter (default 0.6)
        gamma_min (float): Minimum gamma parameter (default 1.0)
        dump_trajectories (bool): Save diffusion trajectories (default False)
        prevalidate_inputs (bool): Validate inputs before design (default False)
        low_memory_mode (bool): Enable low memory consumption mode (default False)
    """
    # Ensure input PDB has valid chain identifiers (RF3 requirement)
    rf3_input_pdb = input_pdb
    if str(input_pdb).endswith('.pdb') and os.path.isfile(input_pdb):
        fixed_pdb_path = f"{output_path}/rf3_input_chainfixed.pdb"
        fill_chain_id = (str(chain_id).strip() if chain_id else "A")[:1]
        replaced_chain_count = 0

        with open(input_pdb, 'r') as fin, open(fixed_pdb_path, 'w') as fout:
            for line in fin:
                if line.startswith(('ATOM', 'HETATM', 'ANISOU', 'TER')) and len(line) >= 22 and line[21] == ' ':
                    line = line[:21] + fill_chain_id + line[22:]
                    replaced_chain_count += 1
                fout.write(line)

        if replaced_chain_count > 0:
            rf3_input_pdb = fixed_pdb_path
            print(f"RF3 input sanitized: set chain ID '{fill_chain_id}' for {replaced_chain_count} records")

    # Create RF3 YAML configuration file
    yaml_path = f"{output_path}/rf3_config.yaml"
    generate_rf3_yaml(
        output_yaml_path=yaml_path,
        input_pdb=rf3_input_pdb,
        contig_map=contig_map,
        ligand_resname=ligand_resname,
        num_designs_batch=num_designs,
        chain_id=chain_id,
        checkpoint_path=checkpoint_path,
        batch_size=batch_size,
        use_classifier_free_guidance=use_classifier_free_guidance,
        cfg_scale=cfg_scale,
        num_timesteps=num_timesteps,
        step_scale=step_scale,
        noise_scale=noise_scale,
        gamma_0=gamma_0,
        gamma_min=gamma_min,
        dump_trajectories=dump_trajectories,
        prevalidate_inputs=prevalidate_inputs,
        low_memory_mode=low_memory_mode
    )
    
    # Build and execute RF3 command
    cmd = (
        f"rfd3 design "
        f"inputs={yaml_path} "
        f"out_dir={output_path} "
        f"ckpt_path={checkpoint_path} "
        f"diffusion_batch_size={batch_size} "
        f"n_batches={num_designs} "
        f"inference_sampler.kind=default "
        f"inference_sampler.cfg_features=[active_donor,active_acceptor,ref_atomwise_rasa] "
        f"inference_sampler.use_classifier_free_guidance={str(use_classifier_free_guidance).lower()} "
        f"inference_sampler.cfg_scale={cfg_scale} "
        f"inference_sampler.center_option=all "
        f"inference_sampler.s_trans=1.0 "
        f"inference_sampler.inference_noise_scaling_factor=1.0 "
        f"inference_sampler.allow_realignment=false "
        f"inference_sampler.num_timesteps={num_timesteps} "
        f"inference_sampler.step_scale={step_scale} "
        f"inference_sampler.noise_scale={noise_scale} "
        f"inference_sampler.p=7 "
        f"inference_sampler.gamma_0={gamma_0} "
        f"inference_sampler.gamma_min={gamma_min} "
        f"inference_sampler.s_jitter_origin=0.0 "
        f"dump_trajectories={str(dump_trajectories).lower()} "
        f"prevalidate_inputs={str(prevalidate_inputs).lower()} "
        f"low_memory_mode={str(low_memory_mode).lower()} "
        f"skip_existing=true"
    )
    cmd_args = cmd.replace("rfd3 ", "", 1)
    
    # Execute RF3 with environment handling
    if path_to_RF3_env:
        env_spec = str(path_to_RF3_env)
        activate_script = os.path.join(env_spec, "bin", "activate")
        rfd3_executable = os.path.join(env_spec, "bin", "rfd3")

        if env_spec.endswith("/bin/activate") and os.path.isfile(env_spec):
            shell(f"source {env_spec} && {cmd}")
        elif os.path.isfile(rfd3_executable):
            shell(f"{rfd3_executable} {cmd_args}")
        elif os.path.isfile(activate_script):
            shell(f"source {activate_script} && {cmd}")
        elif os.path.isdir(env_spec):
            shell(f"conda run -p {env_spec} {cmd}")
        else:
            shell(f"conda run -n {env_spec} {cmd}")
    else:
        shell(cmd)
    
    print(f"RF3 design completed. Output: {output_path}")
