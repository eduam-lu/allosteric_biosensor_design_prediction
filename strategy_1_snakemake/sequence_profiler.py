"""
This is an auxiliary script for the binding site designer pipeline.
Given the reference PDB file, it extracts reference values regarding
- Boltz pLDDT, and affinity
- GNina scores
- DSSP
- Sequence profile regarding: hydrophobicity, SAP score, charge and aliphatic index.

The extracted values are stored in a dictionary that can be used as input for the subsequent steps of the pipeline.
"""
import os
import re
import json
import subprocess
from pathlib import Path
from datetime import datetime

import numpy as np
import yaml
import pandas as pd
import matplotlib.pyplot as plt
import Bio.PDB as PDB
from Bio.PDB import PDBParser, MMCIFParser, DSSP, NeighborSearch, Selection
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.SASA import ShrakeRupley
import biotite.structure.io as bsio

################################################################################
### INPUTS AND PARAMS - EDIT THIS SECTION FOR CONFIGURATION ###
################################################################################

# ============ INPUT/OUTPUT PATHS ============
INPUT_PDB_PATH = "/proj/berzelius-2023-361/users/x_eduam/2P9H.pdb"
OUTPUT_JSON_PATH = "./2P9H_profile/reference_profile.json"

# ============ DIRECTORY PATHS ============
BOLTZ_RUNS_DIR = None                    # Set to directory path if Boltz results exist (e.g., ./boltz_output)
GNINA_OUTPUT_DIR = "./2P9H_profile/gnina_output"

# ============ PROTEIN STRUCTURE PARAMETERS ============
CHAIN_ID = 'A'                            # Default chain ID belonging to the protein
CHAINS_TO_KEEP = ['A']                      # List of chain IDs to keep (e.g., ['A', 'B']). If empty, all chains are kept.
COMPUTE_DSSP = False                     # Compute secondary structure (requires DSSP installed)

# ============ SEQUENCE PROFILE PARAMETERS ============
SEQUENCE_PROFILE_WINDOW_SIZE = 9         # Window size for sliding window profiles (hydrophobicity, charge, aliphatic index)
SAP_SEARCH_RADIUS = 5.0                  # Search radius in Ångströms for SAP score computation

# ============ BOLTZ PARAMETERS ============
BOLTZ_CONFIG = {
    'enabled': True,                    # Set to True to run Boltz sampling
    'ligand_id': 'IPT',                   # 3-letter PDB residue code for ligand (e.g., 'LIG', 'ATP')
    'ligand_smiles': 'CC(C)S[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O',              # SMILES string for ligand
    'n_iterations': 5,                   # Number of Boltz iterations
    'use_msa': True,                     # Enable MSA in Boltz
    'use_forces': False,                  # Use potentials in Boltz
    'no_kernels': True,                 # Disable GPU kernels
    'devices': 1,                        # Number of devices to use
    'recycling_steps': 3,                # Recycling steps
    'sampling_steps': 100,               # Sampling steps for structure prediction
    'diffusion_samples': 1,              # Number of diffusion samples
    'output_format': 'pdb',              # Output format (pdb, mmcif, etc.)
    'sampling_steps_affinity': 100,      # Sampling steps for affinity prediction
    'output_dir': './2P9H_profile/boltz_output',      # Output directory for Boltz runs
    'max_dist': 5.0,                     # Max distance for pocket constraint
    'ligand_search_radius': 5.0,         # Radius to discover binding residues
    'path_to_boltz_env': "/proj/berzelius-2023-361/users/x_eduam/.conda/envs/boltz",           # Path to Boltz conda env (e.g., /home/user/mambaforge/envs/boltz)
}

# ============ GNINA PARAMETERS ============
GNINA_CONFIG = {
    'enabled': True,                    # Set to True to run GNINA scoring
    'gnina_path': '/proj/berzelius-2023-361/users/x_eduam/gnina/build/bin/gnina',      # Path to GNINA executable
    'ligand_id': 'IPT',                  # 3-letter PDB residue code for ligand
    'output_folder': './2P9H_profile/gnina_output',   # Output directory for GNINA results
}

################################################################################
### CONSTANTS - Biochemistry Reference Data ###
################################################################################
HYDROPHOBICITY_SCALE_KYTE_DOOLITTLE = {
    'A': 1.8,
    'C': 2.5,
    'D': -3.5,
    'E': -3.5,
    'F': 2.8,
    'G': -0.4,
    'H': -3.2,
    'I': 4.5,
    'K': -3.9,
    'L': 3.8,
    'M': 1.9,
    'N': -3.5,
    'P': -1.6,
    'Q': -3.5,
    'R': -4.5,
    'S': -0.8,
    'T': -0.7,
    'V': 4.2,
    'W': -0.9,
    'Y': -1.3,
}

CHARGE_SCALE_PH7 = {
    'A': 0.0,
    'C': 0.0,
    'D': -1.0,
    'E': -1.0,
    'F': 0.0,
    'G': 0.0,
    'H': 0.1,
    'I': 0.0,
    'K': 1.0,
    'L': 0.0,
    'M': 0.0,
    'N': 0.0,
    'P': 0.0,
    'Q': 0.0,
    'R': 1.0,
    'S': 0.0,
    'T': 0.0,
    'V': 0.0,
    'W': 0.0,
    'Y': 0.0,
}

AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

HYDROPHOBICITY = {
    'ALA': 0.500, 'ARG': -0.300, 'ASN': -0.500, 'ASP': -0.700, 'CYS': 0.680,
    'GLN': -0.400, 'GLU': -0.600, 'GLY': 0.000, 'HIS': -0.100, 'ILE': 0.943,
    'LEU': 0.943, 'LYS': -0.200, 'MET': 0.738, 'PHE': 1.000, 'PRO': 0.000,
    'SER': 0.300, 'THR': 0.450, 'TRP': 0.849, 'TYR': 0.880, 'VAL': 0.825
}

MAX_SASA_SIDECHAIN = {
    'ALA': 67.0,  'ARG': 196.0, 'ASN': 113.0, 'ASP': 106.0,
    'CYS': 104.0, 'GLN': 144.0, 'GLU': 138.0, 'GLY': 0.1,
    'HIS': 151.0, 'ILE': 140.0, 'LEU': 137.0, 'LYS': 167.0,
    'MET': 160.0, 'PHE': 175.0, 'PRO': 105.0, 'SER': 80.0,
    'THR': 102.0, 'TRP': 217.0, 'TYR': 187.0, 'VAL': 117.0
}

SAP_RADIUS_DEFAULT = 5.0

### FUNCTIONS ###
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



def filter_chains(structure, chains_to_keep=None):
    """
    Filter chains from a structure based on the specified list.
    
    Args:
        structure: Biopython structure object
        chains_to_keep (list, optional): List of chain IDs to keep (e.g., ['A', 'B']). 
                                         If None or empty, returns all chains.
    
    Returns:
        list: List of chain objects that match the filter criteria.
    """
    if not chains_to_keep:
        # Return all chains
        all_chains = []
        for model in structure:
            for chain in model:
                all_chains.append(chain)
        return all_chains
    else:
        # Return only specified chains
        filtered_chains = []
        chains_to_keep_set = set(c.upper() if isinstance(c, str) else c for c in chains_to_keep)
        for model in structure:
            for chain in model:
                if chain.id.upper() in chains_to_keep_set or chain.id in chains_to_keep_set:
                    filtered_chains.append(chain)
        return filtered_chains



def extract_sequence(pdb_path, chain_id='A', chains_to_keep=None):
    """
    Main wrapper to parse file (PDB or CIF) and combine info.
    
    Args:
        pdb_path (str): Path to PDB/CIF file
        chain_id (str): Default chain ID to use if chains_to_keep is None/empty
        chains_to_keep (list, optional): List of chain IDs to include. If provided and non-empty, these chains
                                         are concatenated. If None or empty, uses chain_id.
    
    Returns:
        dict: Sequence information including coordinate_sequence, seqres, etc.
    """
    # 1. Select the correct parser
    if str(pdb_path).endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
        is_cif = True
    else:
        parser = PDBParser(QUIET=True)
        is_cif = False

    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {"error": f"Structure Parsing Error: {str(e)}"}

    # Determine which chains to process
    if chains_to_keep and len(chains_to_keep) > 0:
        # Use chains_to_keep
        chains_to_process = filter_chains(structure, chains_to_keep)
        if not chains_to_process:
            return {"error": f"No chains found matching: {chains_to_keep}"}
    else:
        # Use single chain_id (legacy behavior)
        try:
            chains_to_process = [structure[0][chain_id]]
        except KeyError:
            return {"error": f"Chain {chain_id} not found in structure"}

    # 2. Get SEQRES sequence
    seqres_list = []
    try:
        for chain in chains_to_process:
            seqres = get_seqres_sequence(structure, chain.id, pdb_path)
            if seqres:
                seqres_list.append(seqres)
    except Exception:
        pass

    seqres = ''.join(seqres_list) if seqres_list else None

    # 3. Get Coordinate Sequence
    coord_seq_list = []
    res_count_total = 0
    try:
        for chain in chains_to_process:
            coord_seq, res_count = get_coordinate_sequence(chain)
            coord_seq_list.append(coord_seq)
            res_count_total += res_count
    except Exception as e:
        return {"error": f"Error extracting coordinate sequence: {str(e)}"}

    coord_seq = ''.join(coord_seq_list)

    # 4. Fallback Logic
    if not seqres:
        seqres = coord_seq  # Fallback to coordinate sequence if SEQRES not found/empty

    # 5. Alignment / Padding Logic
    # (If coordinate sequence is shorter than SEQRES, pad the end)
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

    chains_processed = [c.id for c in chains_to_process]

    return {
        "pdb_id": os.path.basename(pdb_path).split('.')[0],
        "seqres": seqres,
        "coordinate_sequence": coord_seq,
        "residue_count": len(coord_seq),
        "found_residues": found_residues,
        "missing_residues": missing_residues,
        "missing_positions": missing_positions,
        "start_residue_number": start_residue_number,
        "sequence_length": len(coord_seq),
        "chains_processed": chains_processed
    }


def get_seqres_sequence(structure, chain_id='A', pdb_path=None):
    """
    Retrieves the SEQRES sequence from the PDB header.
    Converts 3-letter amino acid codes to 1-letter codes.
    Returns None if not found.
    """
    # First try Biopython's parsed header
    if 'seqres' in structure.header:
        seqres_dict = structure.header['seqres']
        if chain_id in seqres_dict:
            seqres_data = seqres_dict[chain_id]

            # Biopython stores SEQRES as a list of 3-letter residue codes
            if isinstance(seqres_data, list):
                seqres_list = seqres_data
            else:
                # If it's already a string, convert it to a list
                seqres_list = str(seqres_data).split()

            # Convert 3-letter codes to 1-letter codes
            one_letter_seq = []
            for res_code in seqres_list:
                res_code = res_code.strip()
                one_letter = AA_3TO1.get(res_code, 'X')  # 'X' for unknown residues
                one_letter_seq.append(one_letter)

            return ''.join(one_letter_seq)

    # Fallback: Parse SEQRES lines directly from PDB file
    # This handles cases where Biopython doesn't populate the header correctly
    if pdb_path:
        try:
            seqres_sequences = {}
            with open(pdb_path, 'r') as f:
                for line in f:
                    if line.startswith('SEQRES'):
                        parts = line.split()
                        if len(parts) >= 5:
                            chain = parts[2]
                            residues = parts[4:]  # Residue codes start from index 4
                            if chain not in seqres_sequences:
                                seqres_sequences[chain] = []
                            seqres_sequences[chain].extend(residues)

            if chain_id in seqres_sequences:
                seqres_list = seqres_sequences[chain_id]

                # Convert 3-letter codes to 1-letter codes
                one_letter_seq = []
                for res_code in seqres_list:
                    res_code = res_code.strip()
                    one_letter = AA_3TO1.get(res_code, 'X')  # 'X' for unknown residues
                    one_letter_seq.append(one_letter)

                return ''.join(one_letter_seq)
        except Exception as e:
            pass

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

def hydrophobicity_profile(sequence, window_size=9, scale=HYDROPHOBICITY_SCALE_KYTE_DOOLITTLE):
    """
    Compute a per-residue hydrophobicity profile using a sliding window.
    Uses the provided hydrophobicity scale and returns one averaged score per
    sequence position.

    Args:
        sequence (str): One-letter amino acid sequence.
        window_size (int): Window width used for averaging. If an even value is
            provided, it is rounded up to the next odd integer.
        scale (dict): Mapping of one-letter amino acid codes to hydrophobicity
            values.

    Returns:
        list[float|None]: Hydrophobicity score per residue. Positions with no
            valid amino acids in the window return None.
    """
    if sequence is None:
        return []

    if window_size < 1:
        raise ValueError("window_size must be >= 1")

    if window_size % 2 == 0:
        window_size += 1

    half_window = window_size // 2
    sequence = sequence.upper()

    profile = []
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window = sequence[start:end]
        window_values = [scale[aa] for aa in window if aa in scale]

        if not window_values:
            profile.append(None)
        else:
            profile.append(sum(window_values) / len(window_values))

    return profile


def _get_sidechain_sasa(residue):
    """Compute sidechain SASA for a residue excluding backbone atoms."""
    sasa = 0.0
    for atom in residue:
        if atom.name in ['N', 'C', 'O', 'CA', 'OXT']:
            continue
        sasa += getattr(atom, 'sasa', 0.0)
    return sasa


def sap_score_profile(pdb_path, radius=SAP_RADIUS_DEFAULT):
    """
    Compute the per-residue SAP profile and return the positive SAP score
    profile that correlates with Rosetta SAP.
    Uses the same solvent accessibility and hydrophobicity logic as
    calculate_sap_score.py.

    Args:
        pdb_path (str): Path to the input PDB file.
        radius (float): Neighbor search radius in Angstroms.

    Returns:
        dict: {
            'pdb_id': str,
            'sequence': str,
            'sap_positive_profile': list[float],
            'residue_count': int
        }
    """
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_path)
    except Exception as e:
        return {
            'pdb_id': os.path.basename(pdb_path).split('.')[0],
            'error': str(e),
            'sequence': '',
            'sap_scores': [],
            'positive_sap_scores': [],
            'residue_count': 0
        }

    residues = []
    atoms_for_search = []

    for model in structure:
        for chain in model:
            for residue in chain:
                if not PDB.is_aa(residue):
                    continue
                if 'CA' not in residue:
                    continue
                residues.append(residue)
                atoms_for_search.append(residue['CA'])

    if not residues:
        return {
            'pdb_id': os.path.basename(pdb_path).split('.')[0],
            'sequence': '',
            'sap_scores': [],
            'positive_sap_scores': [],
            'residue_count': 0
        }

    sr = ShrakeRupley()
    sr.compute(structure, level='A')
    ns = NeighborSearch(atoms_for_search)

    sequence_list = []
    positive_sap_scores = []

    for residue in residues:
        three_letter = residue.get_resname()
        one_letter = AA_3TO1.get(three_letter, 'X')
        sequence_list.append(one_letter)

        target_sap = 0.0
        neighbors = ns.search(residue['CA'].get_coord(), radius, level='R')
        for neighbor in neighbors:
            n_name = neighbor.get_resname()
            actual_area = _get_sidechain_sasa(neighbor)
            max_area = MAX_SASA_SIDECHAIN.get(n_name, 100.0)
            perc_exposed = 0.0
            if max_area >= 1.0:
                perc_exposed = min(actual_area / max_area, 1.0)
            hydrophobicity = HYDROPHOBICITY.get(n_name, 0.0)
            target_sap += perc_exposed * hydrophobicity

        positive_sap_scores.append(target_sap if target_sap > 0 else 0.0)

    return {
        'pdb_id': os.path.basename(pdb_path).split('.')[0],
        'sequence': ''.join(sequence_list),
        'sap_positive_profile': positive_sap_scores,
        'residue_count': len(positive_sap_scores)
    }


def charge_profile(sequence, window_size=9, scale=CHARGE_SCALE_PH7):
    """
    Compute a per-residue charge profile using a sliding window.
    Uses the provided charge scale and returns one averaged charge per
    sequence position.

    Args:
        sequence (str): One-letter amino acid sequence.
        window_size (int): Window width used for averaging. If an even value is
            provided, it is rounded up to the next odd integer.
        scale (dict): Mapping of one-letter amino acid codes to charge values.

    Returns:
        list[float|None]: Average charge per residue. Positions with no valid
            amino acids in the window return None.
    """
    if sequence is None:
        return []

    if window_size < 1:
        raise ValueError("window_size must be >= 1")

    if window_size % 2 == 0:
        window_size += 1

    half_window = window_size // 2
    sequence = sequence.upper()

    profile = []
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window = sequence[start:end]
        window_values = [scale[aa] for aa in window if aa in scale]

        if not window_values:
            profile.append(None)
        else:
            profile.append(sum(window_values) / len(window_values))

    return profile


def aliphatic_index_profile(sequence, window_size=9):
    """
    Compute a per-residue aliphatic index profile using a sliding window.
    Applies AI = 100 * (X_Ala + 2.9 * X_Val + 3.9 * (X_Ile + X_Leu))
    to the residue composition within each window.

    Args:
        sequence (str): One-letter amino acid sequence.
        window_size (int): Window width used for averaging. If an even value is
            provided, it is rounded up to the next odd integer.

    Returns:
        list[float|None]: Aliphatic index per residue. Positions with no valid
            amino acids in the window return None.
    """
    if sequence is None:
        return []

    if window_size < 1:
        raise ValueError("window_size must be >= 1")

    if window_size % 2 == 0:
        window_size += 1

    half_window = window_size // 2
    sequence = sequence.upper()

    profile = []
    for i in range(len(sequence)):
        start = max(0, i - half_window)
        end = min(len(sequence), i + half_window + 1)
        window = sequence[start:end]

        if not window:
            profile.append(None)
            continue

        count_A = window.count('A')
        count_V = window.count('V')
        count_I = window.count('I')
        count_L = window.count('L')
        valid_count = sum(1 for aa in window if aa.isalpha())

        if valid_count == 0:
            profile.append(None)
            continue

        fraction_A = count_A / valid_count
        fraction_V = count_V / valid_count
        fraction_I = count_I / valid_count
        fraction_L = count_L / valid_count

        ai_value = 100.0 * (fraction_A + 2.9 * fraction_V + 3.9 * (fraction_I + fraction_L))
        profile.append(ai_value)

    return profile


def DSSP_profile(input_file):
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

def sequence_profiler(sequence, pdb_path=None, window_size=9, sap_radius=SAP_RADIUS_DEFAULT):
    """
    Compute a dictionary of sequence-derived profiles.

    Args:
        sequence (str): One-letter amino acid sequence.
        pdb_path (str, optional): Path to the PDB file used for SAP profiling.
        window_size (int): Window width used for sliding profiles.
        sap_radius (float): Radius for the SAP neighborhood search.

    Returns:
        dict: {
            'hydrophobicity': list[float|None],
            'charge': list[float|None],
            'aliphatic_index': list[float|None],
            'sap_positive_profile': list[float] | None,
        }
    """
    profiles = {
        'hydrophobicity': hydrophobicity_profile(sequence, window_size=window_size),
        'charge': charge_profile(sequence, window_size=window_size),
        'aliphatic_index': aliphatic_index_profile(sequence, window_size=window_size),
        'sap_positive_profile': None,
    }

    if pdb_path is not None:
        sap_result = sap_score_profile(pdb_path, radius=sap_radius)
        profiles['sap_positive_profile'] = sap_result.get('sap_positive_profile', [])

    return profiles

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

def boltz_sampling(row, ligand_smiles, output_path, pocket_list, max_dist=5.0, n_iterations=1,
                   use_msa=True, use_forces=True, no_kernels=False,
                   path_to_boltz_env=None, devices=1, recycling_steps=3,
                   sampling_steps=100, diffusion_samples=1, output_format='pdb', sampling_steps_affinity=100):
    """
    Run Boltz sampling for a given input row and return the list of iteration folders.

    Args:
        row (dict/Series): Input row containing at least 'file_ID' and optionally sequence.
        ligand_smiles (str): SMILES string for the ligand.
        output_path (str): Root output directory for Boltz runs.
        pocket_list (list): Contact list for the pocket constraint.
        max_dist (float): Pocket max distance.
        n_iterations (int): Number of Boltz runs to execute.
        use_msa (bool): Whether to enable MSA in Boltz.
        use_forces (bool): Whether to enable potentials.
        no_kernels (bool): Whether to disable GPU kernels.
        path_to_boltz_env (str|None): Optional Boltz executable path or conda env path.
        devices (int): Number of devices to pass to Boltz.
        recycling_steps (int): Recycling steps for Boltz.
        sampling_steps (int): Sampling steps for Boltz.
        diffusion_samples (int): Diffusion samples for Boltz.
        output_format (str): Output format passed to Boltz.
        sampling_steps_affinity (int): Affinity sampling steps.

    Returns:
        dict: {
            'run_folders': list[str],
            'processes': list[dict]
        }
    """
    root_path = Path(output_path)
    root_path.mkdir(parents=True, exist_ok=True)

    run_folders = []
    processes = []
    for i in range(1, n_iterations + 1):
        iter_dir = root_path / f"boltz_iter_{i}"
        iter_dir.mkdir(parents=True, exist_ok=True)

        boltz_proc = run_Boltz2(
            row=row,
            ligand_smiles=ligand_smiles,
            output_path=str(iter_dir),
            pocket_list=pocket_list,
            max_dist=max_dist,
            use_msa=use_msa,
            use_forces=use_forces,
            no_kernels=no_kernels,
            path_to_boltz_env=path_to_boltz_env,
            devices=devices,
            recycling_steps=recycling_steps,
            sampling_steps=sampling_steps,
            diffusion_samples=diffusion_samples,
            output_format=output_format,
            sampling_steps_affinity=sampling_steps_affinity
        )

        run_folders.append(str(iter_dir))
        processes.append({
            'iteration': i,
            'output_dir': str(iter_dir),
            'returncode': getattr(boltz_proc, 'returncode', None),
            'stdout': getattr(boltz_proc, 'stdout', None),
            'stderr': getattr(boltz_proc, 'stderr', None)
        })

    return {
        'run_folders': run_folders,
        'processes': processes
    }


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


def run_Boltz2(row,ligand_smiles, output_path,pocket_list, max_dist, use_msa=True, use_forces=True, no_kernels=False,
                path_to_boltz_env="/home/eduardo/mambaforge/envs/boltz",devices=1, recycling_steps=3, sampling_steps=100, diffusion_samples=1, output_format='pdb', sampling_steps_affinity=100):
    # Generate the yaml for this row
    yaml_path = boltz_yaml_generator_w_msa(row, output_path, ligand_smiles, pocket_list, max_dist)

    if path_to_boltz_env:
        boltz_executable = Path(path_to_boltz_env) / "bin" / "boltz"
        boltz_cmd = f"{boltz_executable} predict" if boltz_executable.exists() else "conda run -n boltz boltz predict"
    else:
        boltz_cmd = "conda run -n boltz boltz predict"

    # Base command with mandatory flags
    command_parts = [
        boltz_cmd,
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

def clean_pdb(pdb_path, output_path, keep_ligand_id=None, remove_water=True):
    """
    Clean a PDB file by removing waters and non-target molecules.
    
    Args:
        pdb_path (str): Path to input PDB file.
        output_path (str): Path to save cleaned PDB.
        keep_ligand_id (str, optional): 3-letter residue code to keep (e.g., 'LIG', 'ATP').
                                         If None, keeps only protein.
        remove_water (bool): Whether to remove water molecules. Default: True.
    
    Returns:
        dict: {
            'success': bool,
            'output_path': str,
            'original_residues': int,
            'cleaned_residues': int,
            'removed_count': int,
            'error': str | None
        }
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {
            'success': False,
            'output_path': None,
            'original_residues': 0,
            'cleaned_residues': 0,
            'removed_count': 0,
            'error': f"Error parsing PDB: {str(e)}"
        }
    
    original_count = 0
    cleaned_count = 0
    removed_count = 0
    
    # Build a new structure with only desired residues
    for model in structure:
        residues_to_remove = []
        
        for chain in model:
            for residue in list(chain):
                original_count += 1
                res_name = residue.get_resname().strip()
                
                # Keep protein residues
                if PDB.is_aa(residue):
                    cleaned_count += 1
                    continue
                
                # Keep target ligand if specified
                if keep_ligand_id and res_name == keep_ligand_id:
                    cleaned_count += 1
                    continue
                
                # Remove water if flag is set
                if remove_water and res_name in ['HOH', 'WAT', 'H2O']:
                    residues_to_remove.append((chain, residue))
                    removed_count += 1
                    continue
                
                # Remove other heteroatoms
                residues_to_remove.append((chain, residue))
                removed_count += 1
        
        # Remove marked residues
        for chain, residue in residues_to_remove:
            chain.detach_child(residue.id)
    
    # Write cleaned structure
    try:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        io = PDB.PDBIO()
        io.set_structure(structure)
        io.save(str(output_file))
        
        return {
            'success': True,
            'output_path': str(output_file),
            'original_residues': original_count,
            'cleaned_residues': cleaned_count,
            'removed_count': removed_count,
            'error': None
        }
    except Exception as e:
        return {
            'success': False,
            'output_path': None,
            'original_residues': original_count,
            'cleaned_residues': cleaned_count,
            'removed_count': removed_count,
            'error': f"Error writing PDB: {str(e)}"
        }


def extract_ligand_to_sdf(pdb_path, ligand_id, output_sdf_path):
    """
    Extract ONLY the ligand residue from a PDB file and save as SDF, preserving 3D coordinates.
    
    Args:
        pdb_path (str): Path to input PDB file.
        ligand_id (str): 3-letter residue code for the ligand (e.g., 'LIG', 'ATP').
        output_sdf_path (str): Path to save output SDF file.
    
    Returns:
        dict: {
            'success': bool,
            'output_path': str,
            'ligand_atoms': int,
            'error': str | None
        }
    """
    # Step 1: Use Biopython to extract ONLY the ligand residues
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {
            'success': False,
            'output_path': None,
            'ligand_atoms': 0,
            'error': f"Error parsing PDB: {str(e)}"
        }
    
    ligand_atoms = []
    ligand_residues = []
    
    # Extract only residues matching the ligand_id
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname().strip() == ligand_id:
                    ligand_residues.append(residue)
                    ligand_atoms.extend(residue.get_atoms())
    
    if not ligand_atoms:
        return {
            'success': False,
            'output_path': None,
            'ligand_atoms': 0,
            'error': f"Ligand '{ligand_id}' not found in PDB"
        }
    
    # Step 2: Write ligand-only structure to a temporary PDB
    try:
        temp_pdb = Path(output_sdf_path).parent / f"{Path(output_sdf_path).stem}_ligand_temp.pdb"
        new_structure = PDB.Structure.Structure('ligand')
        new_model = PDB.Model.Model(0)
        new_chain = PDB.Chain.Chain('L')
        
        for residue in ligand_residues:
            new_chain.add(residue.copy())
        
        new_model.add(new_chain)
        new_structure.add(new_model)
        
        io = PDB.PDBIO()
        io.set_structure(new_structure)
        io.save(str(temp_pdb))
        
    except Exception as e:
        return {
            'success': False,
            'output_path': None,
            'ligand_atoms': 0,
            'error': f"Error writing temporary PDB: {str(e)}"
        }
    
    # Step 3: Try to convert to SDF using RDKit
    try:
        from rdkit import Chem
        
        mol = Chem.MolFromPDBFile(str(temp_pdb))
        if mol is None:
            raise Exception("RDKit could not parse ligand PDB")
        
        output_file = Path(output_sdf_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        writer = Chem.SDWriter(str(output_file))
        writer.write(mol)
        writer.close()
        
        # Clean up temporary PDB
        temp_pdb.unlink()
        
        return {
            'success': True,
            'output_path': str(output_file),
            'ligand_atoms': len(ligand_atoms),
            'error': None
        }
    
    except Exception as e:
        # Fallback: Output as PDB if SDF conversion fails
        try:
            output_file = Path(output_sdf_path).parent / f"{Path(output_sdf_path).stem}.pdb"
            # Move temp PDB to final output location
            temp_pdb.rename(output_file)
            
            return {
                'success': True,
                'output_path': str(output_file),
                'ligand_atoms': len(ligand_atoms),
                'error': f'RDKit conversion failed; saved as PDB instead of SDF: {str(e)}'
            }
        except Exception as e2:
            return {
                'success': False,
                'output_path': None,
                'ligand_atoms': len(ligand_atoms),
                'error': f"Error converting to SDF or saving PDB: {str(e2)}"
            }


def gnina_simple(receptor_pdb, ligand_file, output_file, gnina_path='/path/to/gnina'):
    """
    Run GNINA with simple minimization (no autoboxing).
    
    Command: gnina -r rec.pdb -l ligs.sdf --minimize -o minimized.sdf.gz
    
    Args:
        receptor_pdb (str): Path to cleaned receptor PDB.
        ligand_file (str): Path to ligand SDF/PDB file.
        output_file (str): Path to save minimized output.
        gnina_path (str): Path to GNINA executable.
    
    Returns:
        dict: {
            'success': bool,
            'output_path': str,
            'cnn_score': float | None,
            'cnn_affinity': float | None,
            'vina_affinity': float | None,
            'vina_affinity_2': float | None,
            'stdout': str,
            'stderr': str,
            'error': str | None
        }
    """
    gnina_command = f"{gnina_path} -r {receptor_pdb} -l {ligand_file} --minimize -o {output_file}"
    
    try:
        gnina_result = subprocess.run(gnina_command, shell=True, capture_output=True, text=True)
        output = gnina_result.stdout
        
        # Parse output for scores
        match_cnn_score = re.search(r"CNNscore:\s*([-\d.]+)", output)
        CNN_score = float(match_cnn_score.group(1)) if match_cnn_score else None
        
        match_cnn_aff = re.search(r"CNNaffinity:\s*([-\d.]+)", output)
        CNN_affinity = float(match_cnn_aff.group(1)) if match_cnn_aff else None
        
        match_aff = re.search(r"Affinity:\s*([-\d.]+)\s+([-\d.]+)", output)
        if match_aff:
            vina_affinity = float(match_aff.group(1))
            vina_affinity_2 = float(match_aff.group(2))
        else:
            vina_affinity = None
            vina_affinity_2 = None
        
        return {
            'success': True,
            'output_path': output_file,
            'cnn_score': CNN_score,
            'cnn_affinity': CNN_affinity,
            'vina_affinity': vina_affinity,
            'vina_affinity_2': vina_affinity_2,
            'stdout': output,
            'stderr': gnina_result.stderr,
            'error': None
        }
    
    except Exception as e:
        return {
            'success': False,
            'output_path': None,
            'cnn_score': None,
            'cnn_affinity': None,
            'vina_affinity': None,
            'vina_affinity_2': None,
            'stdout': '',
            'stderr': str(e),
            'error': f"GNINA execution error: {str(e)}"
        }



def find_residues_near_ligand(pdb_path, ligand_id, radius=5.0):
    """
    Identify all protein residues within a specified radius of the ligand centroid.
    
    Args:
        pdb_path (str): Path to the PDB file.
        ligand_id (str): 3-letter residue code for the ligand (e.g., 'LIG', 'ATP').
        radius (float): Search radius in Angstroms (default 5.0).
    
    Returns:
        dict: {
            'residues': list[int],  # residue indices/numbers within radius
            'ligand_atoms_count': int,
            'ligand_centroid': list[float],  # [x, y, z]
            'error': str | None
        }
    """
    parser = PDBParser(QUIET=True)
    
    try:
        structure = parser.get_structure('struct', pdb_path)
    except Exception as e:
        return {
            'residues': [],
            'ligand_atoms_count': 0,
            'ligand_centroid': None,
            'error': f"Error parsing PDB: {str(e)}"
        }
    
    # Find all ligand atoms (residue matching ligand_id)
    ligand_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_id:
                    ligand_atoms.extend(residue.get_atoms())
    
    if not ligand_atoms:
        return {
            'residues': [],
            'ligand_atoms_count': 0,
            'ligand_centroid': None,
            'error': f"Ligand '{ligand_id}' not found in structure"
        }
    
    # Calculate ligand centroid
    lig_coords = np.array([atom.get_coord() for atom in ligand_atoms])
    ligand_centroid = np.mean(lig_coords, axis=0).tolist()
    
    # Find protein residues within radius
    residues_in_radius = []
    for model in structure:
        for chain in model:
            for residue in chain:
                # Skip ligands and non-standard residues
                if not PDB.is_aa(residue):
                    continue
                
                # Calculate minimum distance from residue CA to ligand centroid
                if 'CA' in residue:
                    ca_coord = residue['CA'].get_coord()
                    distance = np.linalg.norm(ca_coord - np.array(ligand_centroid))
                    
                    if distance <= radius:
                        residue_number = residue.id[1]
                        residues_in_radius.append(residue_number)
    
    return {
        'residues': sorted(residues_in_radius),
        'ligand_atoms_count': len(ligand_atoms),
        'ligand_centroid': ligand_centroid,
        'radius': radius,
        'error': None
    }


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

### MAIN ###
def generate_reference_profile_json(pdb_path=None, 
                                   output_json_path=None,
                                   boltz_runs_dir=None,
                                   gnina_config=None,
                                   compute_dssp=None,
                                   window_size=None,
                                   sap_radius=None,
                                   chain_id=None,
                                   chains_to_keep=None):
    """
    Orchestrate extraction of all reference profiles and metrics into a comprehensive JSON.
    
    This function uses the configuration defined in the INPUTS AND PARAMS section at the
    top of the script by default. Override any parameter to use custom values.
    
    Args:
        pdb_path (str, optional): Path to the reference PDB structure. Defaults to INPUT_PDB_PATH.
        output_json_path (str, optional): Path to save output JSON. Defaults to OUTPUT_JSON_PATH.
        boltz_runs_dir (str, optional): Directory with Boltz output. Defaults to BOLTZ_RUNS_DIR.
        gnina_config (dict, optional): GNINA configuration. Defaults to GNINA_CONFIG if enabled.
        output_json_path (str, optional): Path to save the output JSON file.
        boltz_runs_dir (str, optional): Directory containing Boltz output runs (e.g., containing boltz_iter_* folders).
        gnina_config (dict, optional): Configuration for GNINA scoring with keys:
            {
                'gnina_path': str,                          # Path to GNINA executable
                'ligand_smiles': str,                       # SMILES string for ligand
                'ligand_id': str,                           # 3-letter PDB residue code (e.g., 'LIG', 'ATP')
                'residues': list,                           # Binding residue indices (AUTO-DISCOVERED if omitted)
                'ligand_search_radius': float,              # Radius in Å to find residues near ligand (default 5.0)
                'size': float,                              # Box size for GNINA scoring
                'exhaustiveness': int (optional),           # GNINA exhaustiveness (default 8)
                'cnn': str (optional),                      # CNN model name (default 'dense')
                'output_folder': str (optional)             # GNINA output directory
            }
            Note: If 'residues' is omitted or empty, and 'ligand_id' is provided,
                  binding residues are automatically discovered within ligand_search_radius.
        compute_dssp (bool): Whether to compute DSSP secondary structure.
        window_size (int): Window size for sequence profiles.
        sap_radius (float): Radius for SAP profile computation.
        chain_id (str): Chain ID to extract (default 'A'). Deprecated, use chains_to_keep instead.
        chains_to_keep (list, optional): List of chain IDs to keep (e.g., ['A', 'B']). 
                                        If None or empty, uses chain_id. If provided, these chains 
                                        are concatenated into a single sequence. Defaults to CHAINS_TO_KEEP.
    
    Returns:
        dict: Comprehensive reference profile dictionary with structure:
        {
            'metadata': {...},
            'sequence_info': {...},
            'sequence_profiles': {...},
            'boltz_metrics': {...} or None,
            'gnina_scores': {...} or None,
            'dssp_info': {...} or None,
            'errors': [...]
        }
    """
    # ============ USE DEFAULTS FROM CONFIGURATION IF NOT PROVIDED ============
    if pdb_path is None:
        pdb_path = INPUT_PDB_PATH
    if output_json_path is None:
        output_json_path = OUTPUT_JSON_PATH
    if boltz_runs_dir is None:
        boltz_runs_dir = BOLTZ_RUNS_DIR
    if gnina_config is None and GNINA_CONFIG.get('enabled'):
        gnina_config = GNINA_CONFIG.copy()
    if compute_dssp is None:
        compute_dssp = COMPUTE_DSSP
    if window_size is None:
        window_size = SEQUENCE_PROFILE_WINDOW_SIZE
    if sap_radius is None:
        sap_radius = SAP_SEARCH_RADIUS
    if chain_id is None:
        chain_id = CHAIN_ID
    if chains_to_keep is None:
        chains_to_keep = CHAINS_TO_KEEP
    
    errors = []
    reference_data = {
        'metadata': {
            'pdb_path': str(pdb_path),
            'chain_id': chain_id,
            'chains_to_keep': chains_to_keep if chains_to_keep else "all",
            'extraction_timestamp': datetime.now().isoformat(),
            'window_size': window_size,
            'sap_radius': sap_radius
        },
        'sequence_info': {},
        'sequence_profiles': {},
        'boltz_metrics': None,
        'gnina_scores': None,
        'dssp_info': None,
        'errors': []
    }
    
    pdb_path = Path(pdb_path)
    if not pdb_path.exists():
        errors.append(f"PDB file not found: {pdb_path}")
        reference_data['errors'] = errors
        return reference_data
    
    # ============ SEQUENCE EXTRACTION ============
    try:
        seq_info = extract_sequence(str(pdb_path), chain_id=chain_id, chains_to_keep=chains_to_keep)
        if 'error' in seq_info:
            errors.append(f"Sequence extraction error: {seq_info['error']}")
        else:
            reference_data['sequence_info'] = seq_info
    except Exception as e:
        errors.append(f"Sequence extraction exception: {str(e)}")
    
    # ============ SEQUENCE PROFILES ============
    if reference_data['sequence_info'] and 'coordinate_sequence' in reference_data['sequence_info']:
        sequence = reference_data['sequence_info']['coordinate_sequence']
        
        try:
            # Compute all sequence-based profiles
            profiles = sequence_profiler(
                sequence=sequence,
                pdb_path=str(pdb_path),
                window_size=window_size,
                sap_radius=sap_radius
            )
            
            reference_data['sequence_profiles'] = {
                'hydrophobicity': profiles.get('hydrophobicity', []),
                'charge': profiles.get('charge', []),
                'aliphatic_index': profiles.get('aliphatic_index', []),
                'sap_positive_profile': profiles.get('sap_positive_profile', []),
            }
        except Exception as e:
            errors.append(f"Sequence profile computation error: {str(e)}")
    
    # ============ BOLTZ METRICS ============
    if BOLTZ_CONFIG.get('enabled', False):
        try:
            # Prepare row data for Boltz sampling
            if reference_data['sequence_info'] and 'coordinate_sequence' in reference_data['sequence_info']:
                sequence = reference_data['sequence_info']['seqres']
                pdb_id = reference_data['sequence_info'].get('pdb_id', 'protein')
                
                # Auto-discover binding residues for pocket constraint
                if BOLTZ_CONFIG.get('ligand_id'):
                    residue_search = find_residues_near_ligand(
                        pdb_path=str(pdb_path),
                        ligand_id=BOLTZ_CONFIG['ligand_id'],
                        radius=BOLTZ_CONFIG.get('ligand_search_radius', 5.0)
                    )
                    
                    if residue_search['error']:
                        errors.append(f"Boltz residue discovery error: {residue_search['error']}")
                        pocket_list = []
                    else:
                        # Convert residue indices to [CHAIN, RES_IDX] format for pocket_list
                        pocket_list = [[chain_id, res_idx] for res_idx in residue_search['residues']]
                else:
                    pocket_list = []
                
                # Create row dict for boltz_sampling
                boltz_row = {
                    'file_ID': pdb_id,
                    'sequence': sequence
                }
                
                # Prepare Boltz output directory
                boltz_output_dir = Path(BOLTZ_CONFIG.get('output_dir', './boltz_output'))
                
                # Run Boltz sampling
                boltz_result = boltz_sampling(
                    row=boltz_row,
                    ligand_smiles=BOLTZ_CONFIG.get('ligand_smiles', 'CCO'),
                    output_path=str(boltz_output_dir),
                    pocket_list=pocket_list,
                    max_dist=BOLTZ_CONFIG.get('max_dist', 5.0),
                    n_iterations=BOLTZ_CONFIG.get('n_iterations', 1),
                    use_msa=BOLTZ_CONFIG.get('use_msa', True),
                    use_forces=BOLTZ_CONFIG.get('use_forces', True),
                    no_kernels=BOLTZ_CONFIG.get('no_kernels', False),
                    path_to_boltz_env=BOLTZ_CONFIG.get('path_to_boltz_env'),
                    devices=BOLTZ_CONFIG.get('devices', 1),
                    recycling_steps=BOLTZ_CONFIG.get('recycling_steps', 3),
                    sampling_steps=BOLTZ_CONFIG.get('sampling_steps', 100),
                    diffusion_samples=BOLTZ_CONFIG.get('diffusion_samples', 1),
                    output_format=BOLTZ_CONFIG.get('output_format', 'pdb'),
                    sampling_steps_affinity=BOLTZ_CONFIG.get('sampling_steps_affinity', 100)
                )
                
                # Process Boltz output - aggregate across all iterations
                if boltz_result['run_folders']:
                    root_boltz_dir = Path(boltz_output_dir)
                    boltz_metrics = process_boltz_output(str(root_boltz_dir), pdb_search_path=str(root_boltz_dir))
                    reference_data['boltz_metrics'] = {
                        'avg_confidence_score': boltz_metrics.get('avg_confidence_score'),
                        'avg_affinity_pred_value': boltz_metrics.get('avg_affinity_pred_value'),
                        'avg_affinity_probability_binary': boltz_metrics.get('avg_affinity_probability_binary'),
                        'avg_plddt_per_position': boltz_metrics.get('avg_plddt_per_position', []),
                        'num_confidence_files': boltz_metrics.get('num_confidence_files', 0),
                        'num_affinity_files': boltz_metrics.get('num_affinity_files', 0),
                        'num_pdbs': boltz_metrics.get('num_pdbs', 0),
                        'source_directory': str(root_boltz_dir),
                        'run_folders': boltz_result['run_folders']
                    }
            else:
                errors.append("Boltz enabled but sequence info not available")
        except Exception as e:
            errors.append(f"Boltz sampling execution error: {str(e)}")
    elif boltz_runs_dir:
        # Fallback: if an existing Boltz output directory is specified, process it
        boltz_path = Path(boltz_runs_dir)
        if boltz_path.exists():
            try:
                boltz_metrics = process_boltz_output(str(boltz_path))
                reference_data['boltz_metrics'] = {
                    'avg_confidence_score': boltz_metrics.get('avg_confidence_score'),
                    'avg_affinity_pred_value': boltz_metrics.get('avg_affinity_pred_value'),
                    'avg_affinity_probability_binary': boltz_metrics.get('avg_affinity_probability_binary'),
                    'avg_plddt_per_position': boltz_metrics.get('avg_plddt_per_position', []),
                    'num_confidence_files': boltz_metrics.get('num_confidence_files', 0),
                    'num_affinity_files': boltz_metrics.get('num_affinity_files', 0),
                    'num_pdbs': boltz_metrics.get('num_pdbs', 0),
                    'source_directory': str(boltz_path)
                }
            except Exception as e:
                errors.append(f"Boltz processing error: {str(e)}")
        else:
            errors.append(f"Boltz runs directory not found: {boltz_path}")
    
    # ============ GNINA SCORES ============
    if gnina_config and gnina_config.get('enabled', False):
        try:
            # Check required keys
            required_keys = ['gnina_path', 'ligand_id']
            if not all(k in gnina_config for k in required_keys):
                errors.append(f"GNINA config missing keys. Required: {required_keys}")
            else:
                gnina_output_dir = Path(gnina_config.get('output_folder', './gnina_output'))
                gnina_output_dir.mkdir(parents=True, exist_ok=True)
                
                # Step 1: Clean the receptor PDB (remove water, keep target ligand)
                cleaned_pdb = gnina_output_dir / f"{pdb_path.stem}_cleaned.pdb"
                clean_result = clean_pdb(
                    pdb_path=str(pdb_path),
                    output_path=str(cleaned_pdb),
                    keep_ligand_id=gnina_config['ligand_id'],
                    remove_water=True
                )
                
                if not clean_result['success']:
                    errors.append(f"GNINA PDB cleaning error: {clean_result['error']}")
                else:
                    reference_data['metadata']['gnina_cleaning'] = {
                        'original_residues': clean_result['original_residues'],
                        'cleaned_residues': clean_result['cleaned_residues'],
                        'removed_count': clean_result['removed_count'],
                        'cleaned_pdb': str(cleaned_pdb)
                    }
                    
                    # Step 2: Extract ligand to SDF (or PDB)
                    ligand_file = gnina_output_dir / f"{pdb_path.stem}_ligand.sdf"
                    extract_result = extract_ligand_to_sdf(
                        pdb_path=str(pdb_path),
                        ligand_id=gnina_config['ligand_id'],
                        output_sdf_path=str(ligand_file)
                    )
                    
                    if not extract_result['success']:
                        errors.append(f"GNINA ligand extraction error: {extract_result['error']}")
                    else:
                        # Step 3: Remove ligand from the cleaned PDB (keep only protein)
                        receptor_only_pdb = gnina_output_dir / f"{pdb_path.stem}_receptor_only.pdb"
                        receptor_clean_result = clean_pdb(
                            pdb_path=str(cleaned_pdb),
                            output_path=str(receptor_only_pdb),
                            keep_ligand_id=None,  # Remove all non-protein residues
                            remove_water=True
                        )
                        
                        if not receptor_clean_result['success']:
                            errors.append(f"GNINA receptor-only PDB creation error: {receptor_clean_result['error']}")
                            # Fallback: use cleaned_pdb with ligand still in it
                            receptor_pdb_to_use = str(cleaned_pdb)
                        else:
                            receptor_pdb_to_use = str(receptor_only_pdb)
                            reference_data['metadata']['gnina_receptor_only'] = {
                                'output_path': str(receptor_only_pdb),
                                'residues_count': receptor_clean_result['cleaned_residues']
                            }
                        
                        # Step 4: Run GNINA with simple minimization
                        output_file = gnina_output_dir / f"{pdb_path.stem}_ligand_minimized.sdf.gz"
                        gnina_result = gnina_simple(
                            receptor_pdb=receptor_pdb_to_use,
                            ligand_file=extract_result['output_path'],
                            output_file=str(output_file),
                            gnina_path=gnina_config['gnina_path']
                        )
                        
                        if gnina_result['success']:
                            reference_data['gnina_scores'] = {
                                'cnn_score': gnina_result['cnn_score'],
                                'cnn_affinity': gnina_result['cnn_affinity'],
                                'vina_affinity': gnina_result['vina_affinity'],
                                'vina_affinity_2': gnina_result['vina_affinity_2'],
                                'output_file': str(output_file),
                                'receptor_pdb': receptor_pdb_to_use,
                                'ligand_extracted': extract_result['output_path'],
                                'ligand_atoms': extract_result['ligand_atoms']
                            }
                        else:
                            errors.append(f"GNINA scoring error: {gnina_result['error']}")
        
        except Exception as e:
            errors.append(f"GNINA execution error: {str(e)}")
    
    # ============ DSSP SECONDARY STRUCTURE ============
    if compute_dssp:
        try:
            dssp_codes = DSSP_profile(str(pdb_path))
            reference_data['dssp_info'] = {
                'secondary_structure': dssp_codes,
                'length': len(dssp_codes)
            }
        except Exception as e:
            errors.append(f"DSSP computation error: {str(e)}")
    
    
    # ============ FINALIZE ============
    reference_data['errors'] = errors
    
    # Save to JSON if output path is provided
    if output_json_path:
        try:
            output_file = Path(output_json_path)
            output_file.parent.mkdir(parents=True, exist_ok=True)
            
            # Helper function to format JSON with pretty-printed objects but compact arrays
            def json_with_compact_arrays(obj, indent=2):
                """
                Serialize to JSON with pretty-printed objects but compact arrays (one line per array).
                """
                def _format(obj, level=0):
                    if isinstance(obj, dict):
                        if not obj:
                            return '{}'
                        current_indent = ' ' * (level * indent)
                        next_indent = ' ' * ((level + 1) * indent)
                        
                        items = []
                        for k, v in obj.items():
                            formatted_v = _format(v, level + 1)
                            items.append(f'{next_indent}"{k}": {formatted_v}')
                        
                        return '{\n' + ',\n'.join(items) + '\n' + current_indent + '}'
                    
                    elif isinstance(obj, list):
                        # Compact format for arrays - all on one line
                        return json.dumps(obj, default=str)
                    
                    else:
                        return json.dumps(obj, default=str)
                
                return _format(obj, level=0)
            
            # Format and save the reference data
            formatted_json = json_with_compact_arrays(reference_data)
            with open(output_file, 'w') as f:
                f.write(formatted_json)
            
            print(f"Reference profile JSON saved to: {output_file}")
        except Exception as e:
            errors.append(f"Error saving JSON to {output_json_path}: {str(e)}")
            reference_data['errors'] = errors
    
    return reference_data


def top_sequence_profiler(df, output_dir, window_size=None, sap_radius=None, 
                          seq_col='sequence', id_col='id', pdb_col=None,
                          compute_sap=False, reference_hydrophobicity=None,
                          reference_charge=None, reference_aliphatic_index=None,
                          reference_sap_profile=None):
    """
    Batch process sequences from a pandas dataframe and export profile dictionaries
    to JSON files. Creates one JSON file per profile type with entries keyed by sequence ID.
    
    Args:
        df (pd.DataFrame): Input dataframe containing sequences.
        output_dir (str): Directory to save output JSON files.
        window_size (int): Window size for sliding window profiles. Defaults to SEQUENCE_PROFILE_WINDOW_SIZE.
        sap_radius (float): Radius for SAP profile computation. Defaults to SAP_SEARCH_RADIUS.
        seq_col (str): Column name containing sequences. Default: 'sequence'.
        id_col (str): Column name containing sequence identifiers. Default: 'id'.
        pdb_col (str, optional): Column name containing PDB file paths. Required if compute_sap=True.
        compute_sap (bool): Whether to compute SAP profiles. Requires pdb_col. Default: False.
        reference_hydrophobicity (list, optional): Reference hydrophobicity profile. If provided, 
            will be added to output with key 'original'.
        reference_charge (list, optional): Reference charge profile. If provided, 
            will be added to output with key 'original'.
        reference_aliphatic_index (list, optional): Reference aliphatic index profile. If provided, 
            will be added to output with key 'original'.
        reference_sap_profile (list, optional): Reference SAP positive profile. If provided, 
            will be added to output with key 'original'.
    
    Returns:
        dict: {
            'profiles_saved': list[str],  # List of saved JSON files
            'entries_processed': int,      # Number of dataframe entries processed
            'errors': list[str]            # List of errors encountered
        }
    
    Output Files:
        - hydrophobicity.json: Hydrophobicity profiles keyed by sequence ID (with 'original' key if reference provided)
        - charge.json: Charge profiles keyed by sequence ID (with 'original' key if reference provided)
        - aliphatic_index.json: Aliphatic index profiles keyed by sequence ID (with 'original' key if reference provided)
        - sap_positive_profile.json: SAP positive profiles (if compute_sap=True; with 'original' key if reference provided)
    """
    # Use defaults if not provided
    if window_size is None:
        window_size = SEQUENCE_PROFILE_WINDOW_SIZE
    if sap_radius is None:
        sap_radius = SAP_SEARCH_RADIUS
    
    # Initialize output dictionaries with reference data if provided
    hydrophobicity_dict = {}
    charge_dict = {}
    aliphatic_index_dict = {}
    sap_positive_profile_dict = {}
    
    # Add reference profiles with 'original' key
    if reference_hydrophobicity is not None:
        hydrophobicity_dict['original'] = reference_hydrophobicity
    if reference_charge is not None:
        charge_dict['original'] = reference_charge
    if reference_aliphatic_index is not None:
        aliphatic_index_dict['original'] = reference_aliphatic_index
    if reference_sap_profile is not None:
        sap_positive_profile_dict['original'] = reference_sap_profile
    
    errors = []
    profiles_saved = []
    
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Validate dataframe columns
    if seq_col not in df.columns:
        errors.append(f"Sequence column '{seq_col}' not found in dataframe")
        return {
            'profiles_saved': [],
            'entries_processed': 0,
            'errors': errors
        }
    
    if id_col not in df.columns:
        errors.append(f"ID column '{id_col}' not found in dataframe")
        return {
            'profiles_saved': [],
            'entries_processed': 0,
            'errors': errors
        }
    
    if compute_sap and pdb_col and pdb_col not in df.columns:
        errors.append(f"PDB column '{pdb_col}' not found in dataframe (required for SAP computation)")
        compute_sap = False
    
    # Process each row in the dataframe
    for idx, row in df.iterrows():
        try:
            seq_id = row[id_col]
            sequence = row[seq_col]
            
            # Skip if sequence is None or empty
            if not sequence or pd.isna(sequence):
                errors.append(f"Row {idx} (ID: {seq_id}): Sequence is empty or None")
                continue
            
            # Convert to string if needed
            sequence = str(sequence)
            
            # Compute hydrophobicity profile
            try:
                hydro_profile = hydrophobicity_profile(sequence, window_size=window_size)
                hydrophobicity_dict[str(seq_id)] = hydro_profile
            except Exception as e:
                errors.append(f"Row {idx} (ID: {seq_id}): Hydrophobicity computation failed: {str(e)}")
            
            # Compute charge profile
            try:
                charge_prof = charge_profile(sequence, window_size=window_size)
                charge_dict[str(seq_id)] = charge_prof
            except Exception as e:
                errors.append(f"Row {idx} (ID: {seq_id}): Charge computation failed: {str(e)}")
            
            # Compute aliphatic index profile
            try:
                ali_profile = aliphatic_index_profile(sequence, window_size=window_size)
                aliphatic_index_dict[str(seq_id)] = ali_profile
            except Exception as e:
                errors.append(f"Row {idx} (ID: {seq_id}): Aliphatic index computation failed: {str(e)}")
            
            # Compute SAP profile if requested and PDB path is available
            if compute_sap and pdb_col:
                try:
                    pdb_path = row[pdb_col]
                    if not pd.isna(pdb_path) and Path(pdb_path).exists():
                        sap_result = sap_score_profile(pdb_path, radius=sap_radius)
                        if 'error' not in sap_result:
                            sap_positive_profile_dict[str(seq_id)] = sap_result.get('sap_positive_profile', [])
                        else:
                            errors.append(f"Row {idx} (ID: {seq_id}): SAP computation error: {sap_result['error']}")
                    else:
                        errors.append(f"Row {idx} (ID: {seq_id}): PDB file not found or invalid path")
                except Exception as e:
                    errors.append(f"Row {idx} (ID: {seq_id}): SAP computation exception: {str(e)}")
        
        except Exception as e:
            errors.append(f"Row {idx}: Unexpected error during processing: {str(e)}")
            continue
    
    # Save dictionaries to JSON files
    try:
        hydrophobicity_file = output_path / "hydrophobicity.json"
        with open(hydrophobicity_file, 'w') as f:
            json.dump(hydrophobicity_dict, f, indent=2, default=str)
        profiles_saved.append(str(hydrophobicity_file))
        print(f"✓ Saved hydrophobicity profiles: {hydrophobicity_file} ({len(hydrophobicity_dict)} entries)")
    except Exception as e:
        errors.append(f"Error saving hydrophobicity.json: {str(e)}")
    
    try:
        charge_file = output_path / "charge.json"
        with open(charge_file, 'w') as f:
            json.dump(charge_dict, f, indent=2, default=str)
        profiles_saved.append(str(charge_file))
        print(f"✓ Saved charge profiles: {charge_file} ({len(charge_dict)} entries)")
    except Exception as e:
        errors.append(f"Error saving charge.json: {str(e)}")
    
    try:
        aliphatic_file = output_path / "aliphatic_index.json"
        with open(aliphatic_file, 'w') as f:
            json.dump(aliphatic_index_dict, f, indent=2, default=str)
        profiles_saved.append(str(aliphatic_file))
        print(f"✓ Saved aliphatic index profiles: {aliphatic_file} ({len(aliphatic_index_dict)} entries)")
    except Exception as e:
        errors.append(f"Error saving aliphatic_index.json: {str(e)}")
    
    if compute_sap and sap_positive_profile_dict:
        try:
            sap_file = output_path / "sap_positive_profile.json"
            with open(sap_file, 'w') as f:
                json.dump(sap_positive_profile_dict, f, indent=2, default=str)
            profiles_saved.append(str(sap_file))
            print(f"✓ Saved SAP positive profiles: {sap_file} ({len(sap_positive_profile_dict)} entries)")
        except Exception as e:
            errors.append(f"Error saving sap_positive_profile.json: {str(e)}")
    
    return {
        'profiles_saved': profiles_saved,
        'entries_processed': len(hydrophobicity_dict),
        'errors': errors
    }


def plot_profile(profile_dict, xlabel, ylabel, title, output):
    """
    Plot all entries from a profile dictionary as lines on a single plot.
    
    The 'original' entry (if present) is plotted with a thicker linewidth for emphasis.
    All other entries are plotted with normal linewidth.
    
    Args:
        profile_dict (dict): Dictionary with profile identifiers as keys and profile lists/arrays as values.
            Example: {'original': [1.5, 2.0, 1.8, ...], 'seq_1': [1.2, 1.9, 1.7, ...], ...}
        xlabel (str): Label for the x-axis.
        ylabel (str): Label for the y-axis.
        title (str): Title for the plot.
        output (str): Path to save the output plot image file.
    
    Returns:
        dict: {
            'plot_saved': str,   # Path to saved plot
            'entries_plotted': int,  # Number of entries plotted
            'errors': list[str]   # Any errors encountered
        }
    """
    errors = []
    entries_plotted = 0
    
    try:
        # Create figure and axis
        fig, ax = plt.subplots(figsize=(12, 7))
        
        # Plot each entry in the dictionary
        for key, profile in profile_dict.items():
            if profile is None or (isinstance(profile, list) and len(profile) == 0):
                errors.append(f"Skipping entry '{key}': empty or None profile")
                continue
            
            # Convert to numpy array for consistent handling
            profile_array = np.array(profile, dtype=float)
            
            # Skip if array is empty or all NaN
            if profile_array.size == 0 or np.all(np.isnan(profile_array)):
                errors.append(f"Skipping entry '{key}': contains no valid numeric data")
                continue
            
            # Determine linewidth: thicker for 'original', normal for others
            linewidth = 3.0 if key == 'original' else 1.5
            
            # Plot the profile
            x_positions = np.arange(len(profile_array))
            ax.plot(x_positions, profile_array, label=key, linewidth=linewidth, alpha=0.8)
            
            entries_plotted += 1
        
        # Configure plot
        ax.set_xlabel(xlabel, fontsize=12, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
        ax.grid(True, alpha=0.3)
        
        # Tight layout to prevent label cutoff
        plt.tight_layout()
        
        # Save the plot
        output_path = Path(output)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"✓ Plot saved to: {output_path} ({entries_plotted} entries plotted)")
        
        return {
            'plot_saved': str(output_path),
            'entries_plotted': entries_plotted,
            'errors': errors
        }
    
    except Exception as e:
        errors.append(f"Error creating plot: {str(e)}")
        plt.close()
        return {
            'plot_saved': None,
            'entries_plotted': 0,
            'errors': errors
        }


################################################################################
### ENTRY POINT ###
################################################################################

def main():
    """
    Main entry point. Runs the reference profile extraction using configuration
    from the INPUTS AND PARAMS section at the top of the script.
    
    Call this function to execute the extraction pipeline with all configured parameters.
    """
    print("=" * 80)
    print("REFERENCE PROFILE EXTRACTION - BINDING SITE DESIGNER PIPELINE")
    print("=" * 80)
    print()
    
    print("Configuration Summary:")
    print(f"  Input PDB:          {INPUT_PDB_PATH}")
    print(f"  Output JSON:        {OUTPUT_JSON_PATH}")
    print(f"  Chain ID:           {CHAIN_ID}")
    print(f"  Window Size:        {SEQUENCE_PROFILE_WINDOW_SIZE}")
    print(f"  SAP Radius:         {SAP_SEARCH_RADIUS} Å")
    print(f"  Compute DSSP:       {COMPUTE_DSSP}")
    print(f"  Boltz Enabled:      {BOLTZ_CONFIG.get('enabled', False)}")
    print(f"  GNINA Enabled:      {GNINA_CONFIG.get('enabled', False)}")
    print()
    
    # Run the extraction
    try:
        reference_data = generate_reference_profile_json()
        
        print("=" * 80)
        print("EXTRACTION COMPLETE")
        print("=" * 80)
        print()
        
        # Print summary
        if reference_data['errors']:
            print(f"Warnings/Errors ({len(reference_data['errors'])}):")
            for err in reference_data['errors']:
                print(f"  ⚠ {err}")
            print()
        
        print("Extracted Data Summary:")
        if reference_data['sequence_info']:
            seq = reference_data['sequence_info']
            print(f"  Sequence Length:    {seq.get('sequence_length', 'N/A')}")
            print(f"  Found Residues:     {seq.get('found_residues', 'N/A')}")
            print(f"  Missing Positions:  {len(seq.get('missing_positions', []))} gaps")
        
        if reference_data['boltz_metrics']:
            boltz = reference_data['boltz_metrics']
            print(f"  Boltz Confidence:   {boltz.get('avg_confidence_score', 'N/A'):.3f}" if boltz.get('avg_confidence_score') else "  Boltz Confidence:   N/A")
            print(f"  Boltz PDBs Found:   {boltz.get('num_pdbs', 0)}")
        
        if reference_data['gnina_scores']:
            gnina = reference_data['gnina_scores']
            print(f"  GNINA CNN Score:    {gnina.get('cnn_score', 'N/A')}")
            print(f"  Binding Residues:   {len(gnina.get('residues_used', []))} residues")
        
        print()
        print(f"✓ Profile saved to: {OUTPUT_JSON_PATH}")
        
        return reference_data
        
    except Exception as e:
        print(f"✗ ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        return None


if __name__ == '__main__':
    main()

