"""
binding_site_designer.py --structure --ligand --output --help

DESCRIPTION

This script, given an input structure with a ligand bound to it, performs a binding site redesign around a given new ligand, in order to accept it.
Binding site redesign is performed through RF diffusion all atom and LigandMPNN. Candidates are filtered based on prediction, 3D and binding metrics

FLAGS

--structure:

--ligand:

--output:


--help:

WORKFLOW

NOTES
- 

REQUIREMENTS
- Python 3.8+
- Packages: numpy, pandas, biopython, rdkit, pymol, logging
- External tools: Rosetta, RF diffusion, LigandMPNN
- Mamba and environment with SE3nv installed for RF diffusion execution
Eduardo Amo González
2025-2026
"""
### IMPORT MODULES ###########################################################################################################################

from ligand_sampling import run_ligand_sampling_pipeline
import functions_bsd as func
import argparse
import pandas
import json
from pathlib import Path
import sys
import os
import logging

### PARAMS ####################################################################################################################################

### Ligand sampling params
ligand_smiles = "CC1=CC=CC=C1"  # SMILES string of the ligand to be sampled
ligand_name = "Tc"  # Name of the ligand (used in PDB files)

num_conformers = 2  # Number of conformers to generate for the ligand
num_positions = 2  # Number of different positions to sample for the ligand around the protein
conformer_rmsd_cutoff = 0.75  # RMSD cutoff for filtering conformers. Anything below this value will be considered redundant.

rotation = 15  # Rotation angle in degrees for exploring ligand placement
translation = 0.3  # Translation distance in angstroms for exploring ligand placement

# Active site definition params
first_shell_distance = 5.0  # Distance in angstroms to define the first shell of residues around the ligand
second_shell_distance = 5.0  # Distance in angstroms to define the second shell of residues around the first shell
include_second_shell = False  # If true, residues in the second shell will also be included in the redesign
user_defined_active_site = None  # If provided, this list of residues will be used as the active site instead of calculating one
user_defined_residues = None  # If provided, the program ensures this list of residues are included in the active site

### Redesign conditions params

ligand_MPNN_only = False  # If true, only the ligandMPNN step will be performed, skipping RF diffusion
RF_model = "all_atom"  # RF diffusion model to use: "all_atom", "RF1", "RF3"

# RF redesign conditions
consertative_redesign = False  # If true,residues within a secondary structure element will not be redesigned. Except if they are the first of final residue of the element.
segment_extension = 1  # Number of residues to extend the redesign segment on each side
n_termini_extension = 1  # Number of residues to extend the redesign segment at the N-terminus
c_termini_extension = None  # Number of residues to extend the redesign segment at the C-terminus. If None, it will extend to the end of the chain.
user_defined_contig_map = None  # If provided, this contig map will be used for RF diffusion redesign instead of calculating one based on the structure


### RF diffusion execution params
path_to_RFAA_apptainer ="/home/eduardo/allostery/rf_diffusion_all_atom/rf_se3_diffusion.sif"  # Path to the RF diffusion all atom apptainer image
path_to_RFAA_script ="/home/eduardo/allostery/rf_diffusion_all_atom/run_inference.py"  # Path to the RF diffusion inference script
path_to_RFAA_weights ="/home/eduardo/allostery/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt"  # Path to the RF diffusion all atom weights file
path_to_RF1_script= "/home/eduardo/programs/RFdiffusion/scripts/run_inference.py"  # Path to the RF1 inference script
path_to_RF1_env="/home/sofia/miniforge3/envs/SE3-nvidia" # Path to the environment needed to use RF1


num_designs = 1  # Number of designs to generate with RF diffusion
T = 50  # Number of diffusion steps to use in RF diffusion
RFAA_ligand_name = "UNL"  # Name of the ligand to be used in
design_startnum = 1  # Starting number for design numbering
deterministic = True  # If true, RF diffusion will run in deterministic mode for reproducibility

### Ligand MPNN execution params
path_to_ligand_MPNN_script= "/home/eduardo/LigandMPNN/run.py"  # Path to the Ligand MPNN script 
path_to_ligand_MPNN_env= "/home/eduardo/miniforge3/envs/ligandmpnn_env"  # Path to the environment needed to use Ligand MPNN
mpnn_model_type= "ligand_mpnn"  # Type of Ligand MPNN model to use: "base" or "large" 
path_to_mpnn_model = "/home/eduardo/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"

MPNN_num_designs = 10  # Number of designs to generate with Ligand MPNN per input structure
n_batches = 1  # Number of batches to split the MPNN designs into
mpnn_temperature= 0.05  # Sampling temperature for Ligand MPNN 
bias_aa_global= None  # Global amino acid bias for Ligand MPNN
omit_aa_global= None  # Global amino acids to omit for Ligand MPNN
side_chain_context = 0  # 0 or 1, Ligand MPNN will use side chain context information during design
first_shell_only = False  # If true, only residues in the first shell will be redesigned by Ligand MPNN

user_defined_mpnn_redesign = None  # If provided, this list of residues will be used as the redesign list for Ligand MPNN

top_n_mpnn_candidates= 5

### 1D filtering
filter_1d_window_size = 10
filter_1d_treshold = 10

### ESM
path_to_ESM_env='/home/eduardo/miniforge3/envs/esm_fold_env'
path_to_ESM_script='/home/eduardo/allostery/strategy_1/esm_high_throughput.py'
### INPUT CHECK ###############################################################################################################################

parser = argparse.ArgumentParser(
    description="Binding Site Designer: Redesign binding sites around a new ligand using RF diffusion and LigandMPNN",
)
parser.add_argument('--structure', help="Path to the input structure file (PDB format)", type=str, required=True)
parser.add_argument('--ligand', help="SMILES string of the ligand to be sampled", type=str, required=True)
parser.add_argument('--output', help="Path to the output folder where results will be saved", type=str, required=True)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# Generate output directory if it doesn't exist
Path(args.output).mkdir(parents=True, exist_ok=True)

### SET UP LOG FILE #############################################################################################################################
logging.basicConfig(
    filename=f"{str(args.output)}/binding_site_designer.log",       # Log file path
    filemode='a',                   # Append mode
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)

def save_to_log(message):
    logging.info(message)

### MAIN EXECUTION ##############################################################################################################################

### Generate new pdbs with the new ligand

run_ligand_sampling_pipeline(
        ligand_smiles=args.ligand,
        structure_path=args.structure,
        output_path=f"{args.output}/ligand_sampling",
        num_conformers=num_conformers,
        num_positions=num_positions,
        ligand_name=ligand_name, # Keep default or add CLI arg for this,
        rotation=rotation,
        translation=translation,
        conformer_rmsd_cutoff=conformer_rmsd_cutoff
        )

### Generate redesign information

pdb_info, printable_pdb_info = func.extract_pdb_info(pdb_path=args.structure,first_shell_distance=first_shell_distance,
                                 second_shell_distance=second_shell_distance,
                                 user_defined_active_site=user_defined_active_site,
                                 user_defined_residues=user_defined_residues)
# Save pdb info to json

with open(f"{args.output}/pdb_info.json", "w") as f:
    json.dump(printable_pdb_info, f, indent=2)

# Generate contig map info

contig_map,segments = func.list_to_contig_map(
    chain_id = pdb_info["chain_id"],
    seq_length= pdb_info["sequence_length"],
    start = pdb_info["start_residue_number"],
    missing_residues= pdb_info["missing_positions"],
    active_site = pdb_info["active_site"],
    segment_extension = segment_extension,
    n_termini_extension = n_termini_extension,
    c_termini_extension = c_termini_extension,
    conservative_RF = consertative_redesign,
    DSSP_string = pdb_info["dssp_string"],
    user_defined_contig_map = user_defined_contig_map
    )

print(contig_map)
pdb_info["contig_map"] = contig_map
pdb_info["redesign_segments"] = segments


### RF Diffusion module

if not ligand_MPNN_only:
    if RF_model == "all_atom":
        # Skip if output already exists and is not empty
        if os.path.exists(f"{args.output}/RFallatom_designs") and os.listdir(f"{args.output}/RFallatom_designs"):
            print(f"RFallatom designs already exist in {args.output}/RFallatom_designs. Skipping RFallatom execution.")
        else:
            print("Running RF all atom redesign...")
            for pdb_file in Path(f"{args.output}/ligand_sampling/final_pdbs").glob("*.pdb"):
                func.run_rfAA(
                    output_path=args.output,
                    input_pdb = pdb_file, # Change
                    contig_map=contig_map,
                    num_designs=num_designs,
                    T=T,
                    path_to_RFAA_apptainer = path_to_RFAA_apptainer,
                    path_to_RFAA_script = path_to_RFAA_script ,
                    path_to_RFAA_weights=path_to_RFAA_weights ,
                    inference_ligand=RFAA_ligand_name,
                    design_startnum=design_startnum,
                    deterministic=deterministic
                )
            func.move_pdbs_to_folder(
                input_folder=f"{args.output}/RFallatom_designs",
                output_folder=f"{args.output}/RF_final_pdbs")
    elif RF_model == "RF3":
        if os.path.exists(f"{args.output}/RF3_designs") and os.listdir(f"{args.output}/RF3_designs"):
            print(f"RF3 designs already exist in {args.output}/RF3_designs. Skipping RF3 execution.")
        else:
            print("Running RF3 redesign...")
            for pdb_file in Path(f"{args.output}/ligand_sampling/final_pdbs").glob("*.pdb"):
                func.run_rf3()
            func.move_pdbs_to_folder(
                input_folder=f"{args.output}/RF3_designs",
                output_folder=f"{args.output}/RF_final_pdbs")
    elif RF_model == "RF1":
        print("Running RF1 redesign...")
        # Skip if output already exists and is not empty
        if os.path.exists(f"{args.output}/RF1_designs") and os.listdir(f"{args.output}/RF1_designs"):
            print(f"RF1 designs already exist in {args.output}/RF1_designs. Skipping RF1 execution.")
        else:
            for pdb_file in Path(f"{args.output}/ligand_sampling/final_pdbs").glob("*.pdb"):
                func.run_rf1(
                    output_path=args.output,
                    input_pdb = pdb_file, #Changes
                    contig_map=contig_map,
                    num_designs=num_designs,
                    T=T,
                    path_to_RF1_script= path_to_RF1_script,
                    path_to_RF1_env=path_to_RF1_env
                )
            func.move_pdbs_to_folder(
                input_folder=f"{args.output}/RF1_designs",
                output_folder=f"{args.output}/RF_final_pdbs")
    else:
        print("Invalid RF diffusion model specified. Choose either 'all_atom', 'RF1', or 'RF3'.")
        sys.exit(1)



### Ligand MPNN

# Generate JSONS

if ligand_MPNN_only:
    mpnn_pdb_folder = f"{args.output}/ligand_sampling/final_pdbs"
else:
    mpnn_pdb_folder = f"{args.output}/RF_final_pdbs"

mpnn_output_folder = f"{args.output}/MPNN_designs"
mpnn_json_path = f"{mpnn_output_folder}/MPNN_jsons"
Path(mpnn_json_path).mkdir(parents=True, exist_ok=True)

func.generate_MPNN_jsons(pdb_folder= mpnn_pdb_folder, 
                         output_json = mpnn_json_path,
                         user_redesign_list = user_defined_mpnn_redesign, 
                         first_shell = pdb_info["first_shell"], 
                         chain_id = pdb_info["chain_id"], 
                         original_seq = pdb_info["sequence_seqres"],
                         segments=pdb_info["redesign_segments"],
                         start_position= pdb_info["start_residue_number"],
                         first_shell_only = first_shell_only)

# Run ligand MPNN

func.run_ligand_MPNN(output_folder= mpnn_output_folder,
                    num_designs = MPNN_num_designs, 
                    n_batches = n_batches, 
                    path_to_ligand_MPNN_script = path_to_ligand_MPNN_script, 
                    path_to_ligand_MPNN_env = path_to_ligand_MPNN_env,
                    json_path = mpnn_json_path,
                    model_type = mpnn_model_type, 
                    temperature = mpnn_temperature, 
                    bias_aa_global = bias_aa_global, 
                    omit_aa_global = omit_aa_global, 
                    side_chain_context = side_chain_context,
                    model_path=path_to_mpnn_model)

# Generate MPNN dataframe

MPNN_df = func.process_MPNN_folder(
    folder = f"{args.output}/MPNN_designs/seqs",
    top_n = top_n_mpnn_candidates
)
print(MPNN_df.head(5))

### High throughput filtering of candidates #########################################################

### 1D filter
# Polyalanine filter
MPNN_df = func.filter_dataframe_1D(MPNN_df,window_size=filter_1d_window_size,threshold=filter_1d_treshold, seq_col='seq',aa ='A')
# Polyglutamate filter
MPNN_df = func.filter_dataframe_1D(MPNN_df,window_size=filter_1d_window_size,threshold=filter_1d_treshold, seq_col='seq',aa ='E')

MPNN_df.to_csv(f'{args.output}/MPNN_df.csv')
### High throughput structure prediction with ESM fold
func.run_ESMfold(f'{args.output}/MPNN_df.csv',f'{args.output}/ESM_predictions/',path_to_ESM_env=path_to_ESM_env,path_to_ESM_script=path_to_ESM_script)
### First 3D filter (3D quality)

### Chai prediction 

### Second 3D filter (3D quality and pocket metrics)

### Final filtering and selection of best candidates

