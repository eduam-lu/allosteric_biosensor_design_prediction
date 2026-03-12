"""
binding_site_designer.py --structure --ligand --output --help

DESCRIPTION

This script, given an input structure with a ligand bound to it, performs a binding site redesign around a given new ligand, in order to accept it.
Binding site redesign is performed through RF diffusion all atom and LigandMPNN. Candidates are filtered based on prediction, 3D and binding metrics

FLAGS

--config: Path to the YAML configuration file with all parameters for the redesign process. See below for details on the expected format of this file.

--output: Path to the output folder where results will be saved. The script will create this folder if it doesn't exist.

--help: Show a brief help message and exit.

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
import pandas as pd
import json
from pathlib import Path
import sys
import os
import logging
import yaml

### INPUT CHECK ###############################################################################################################################
parser = argparse.ArgumentParser(
    description="Binding Site Designer: Redesign binding sites around a new ligand using RF diffusion and LigandMPNN",
)
parser.add_argument('--output', help="Path to the output folder where results will be saved", type=str, required=True)
parser.add_argument('--config', help="Path to the YAML configuration file", type=str, required=True)
parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

args = parser.parse_args()

# Generate output directory if it doesn't exist
Path(args.output).mkdir(parents=True, exist_ok=True)

### LOAD PARAMS FROM CONFIG ####################################################################################################################

with open(args.config, 'r') as file:
    config = yaml.safe_load(file)

### Basic params
# - structure: Path to the input structure file (PDB format)
# - ligand_smiles: SMILES string of the ligand to be sampled
# - ligand_name: Name of the ligand (used in PDB files)
structure_path = config.get('structure_path')
ligand_smiles = config.get('ligand_smiles')
ligand_name = config.get('ligand_name', "Tc")

### Ligand sampling params
# - num_conformers: Number of conformers to generate for the ligand
# - num_positions: Number of different positions to sample for the ligand around the protein
# - conformer_rmsd_cutoff: RMSD cutoff for filtering conformers. Anything below this value will be considered redundant
# - rotation: Rotation angle in degrees for exploring ligand placement
# - translation: Translation distance in angstroms for exploring ligand placement
num_conformers = config.get('num_conformers', 2)
num_positions = config.get('num_positions', 1)
conformer_rmsd_cutoff = config.get('conformer_rmsd_cutoff', 0.75)
rotation = config.get('rotation', 15)
translation = config.get('translation', 0.3)

### Active site definition params
# - first_shell_distance: Distance in angstroms to define the first shell of residues around the ligand
# - second_shell_distance: Distance in angstroms to define the second shell of residues around the first shell
# - include_second_shell: If true, residues in the second shell will also be included in the redesign
# - user_defined_active_site: If provided, this list of residues will be used as the active site instead of calculating one
# - user_defined_residues: If provided, the program ensures this list of residues are included in the active site
first_shell_distance = config.get('first_shell_distance', 5.0)
second_shell_distance = config.get('second_shell_distance', 5.0)
include_second_shell = config.get('include_second_shell', False)
user_defined_active_site = config.get('user_defined_active_site', None)
user_defined_residues = config.get('user_defined_residues', None)

### Redesign conditions params
# - ligand_MPNN_only: If true, only the ligandMPNN step will be performed, skipping RF diffusion
# - RF_folder: If provided, this folder with RF diffusion designs will be used as input for LigandMPNN instead of running RF diffusion
# - RF_model: RF diffusion model to use ("all_atom", "RF1", "RF3")
# - consertative_redesign: If true, residues within a secondary structure element will not be redesigned. Except if they are the first or final residue of the element (NOT IMPLEMENTED YET)
# - segment_extension: Number of residues to extend the redesign segment on each side
# - n_termini_extension: Number of residues to extend the redesign segment at the N-terminus
# - c_termini_extension: Number of residues to extend the redesign segment at the C-terminus. If None, it will extend to the end of the chain
# - user_defined_contig_map: If provided, this contig map will be used for RF diffusion redesign instead of calculating one
ligand_MPNN_only = config.get('ligand_MPNN_only', False)
RF_folder = config.get('RF_folder', None)
RF_model = config.get('RF_model', "all_atom")
consertative_redesign = config.get('consertative_redesign', False)
segment_extension = config.get('segment_extension', 1)
n_termini_extension = config.get('n_termini_extension', 1)
c_termini_extension = config.get('c_termini_extension', None)
user_defined_contig_map = config.get('user_defined_contig_map', None)

### RF diffusion execution params
# - path_to_RFAA_apptainer: Path to the RF diffusion all atom apptainer image
# - path_to_RFAA_script: Path to the RF diffusion inference script
# - path_to_RFAA_weights: Path to the RF diffusion all atom weights file
# - path_to_RF1_script: Path to the RF1 inference script
# - path_to_RF1_env: Path to the environment needed to use RF1
# - num_designs: Number of designs to generate with RF diffusion
# - T: Number of diffusion steps to use in RF diffusion
# - RFAA_ligand_name: Name of the ligand to be used in RFAA
# - design_startnum: Starting number for design numbering
# - deterministic: If true, RF diffusion will run in deterministic mode for reproducibility
path_to_RFAA_apptainer = config.get('path_to_RFAA_apptainer')
path_to_RFAA_script = config.get('path_to_RFAA_script')
path_to_RFAA_weights = config.get('path_to_RFAA_weights')
path_to_RF1_script = config.get('path_to_RF1_script')
path_to_RF1_env = config.get('path_to_RF1_env')
num_designs = config.get('num_designs', 1)
T = config.get('T', 50)
RFAA_ligand_name = config.get('RFAA_ligand_name', "UNL")
design_startnum = config.get('design_startnum', 1)
deterministic = config.get('deterministic', True)

### Ligand MPNN execution params
# - path_to_ligand_MPNN_script: Path to the Ligand MPNN script
# - path_to_ligand_MPNN_env: Path to the environment needed to use Ligand MPNN
# - mpnn_model_type: Type of Ligand MPNN model to use ("base" or "ligand_mpnn")
# - path_to_mpnn_model: Path to the ligand MPNN model weights
# - MPNN_num_designs: Number of designs to generate with Ligand MPNN per input structure
# - n_batches: Number of batches to split the MPNN designs into
# - mpnn_temperature: Sampling temperature for Ligand MPNN
# - bias_aa_global: Global amino acid bias for Ligand MPNN
# - omit_aa_global: Global amino acids to omit for Ligand MPNN
# - side_chain_context: 0 or 1, Ligand MPNN will use side chain context information during design
# - first_shell_only: If true, only residues in the first shell will be redesigned by Ligand MPNN
# - user_defined_mpnn_redesign: If provided, this list of residues will be used as the redesign list for Ligand MPNN
# - top_n_mpnn_candidates: Number of top MPNN sequence designs to keep and evaluate
path_to_ligand_MPNN_script = config.get('path_to_ligand_MPNN_script')
path_to_ligand_MPNN_env = config.get('path_to_ligand_MPNN_env')
mpnn_model_type = config.get('mpnn_model_type', "ligand_mpnn")
path_to_mpnn_model = config.get('path_to_mpnn_model')
MPNN_num_designs = config.get('MPNN_num_designs', 10)
n_batches = config.get('n_batches', 1)
mpnn_temperature = config.get('mpnn_temperature', 0.05)
bias_aa_global = config.get('bias_aa_global', None)
omit_aa_global = config.get('omit_aa_global', None)
side_chain_context = config.get('side_chain_context', 0)
first_shell_only = config.get('first_shell_only', False)
user_defined_mpnn_redesign = config.get('user_defined_mpnn_redesign', None)
top_n_mpnn_candidates = config.get('top_n_mpnn_candidates', 5)

### ESM
# - path_to_ESM_env: Path to the ESM fold environment
# - path_to_ESM_script: Path to the high throughput ESM folding python script
# - path_to_ESM_image: Path to the ESM fold apptainer/singularity image
path_to_ESM_env = config.get('path_to_ESM_env')
path_to_ESM_script = config.get('path_to_ESM_script')
path_to_ESM_image = config.get('path_to_ESM_image')

### 1D filtering
# - filter_1d_window_size: Window size for detecting single-amino acid repeats (like Poly-A or Poly-E)
# - filter_1d_treshold: Threshold amount of repeating amino acids to trigger filtering a sequence
filter_1d_window_size = config.get('filter_1d_window_size', 10)
filter_1d_treshold = config.get('filter_1d_treshold', 10)

### First 3D filtering
# - clash_distance: Distance threshold (Angstroms) to consider two atoms clashing
# - bond_distance: Distance threshold (Angstroms) to consider two atoms bonded
# - MIN_PLDDT_1: Minimum average pLDDT score required to pass the first filter
# - MIN_RMSD_1: Minimum RMSD to ensure the design has structurally diverged enough from the original
# - MAX_RMSD_1: Maximum RMSD to ensure the design hasn't lost overall structural integrity
# - MAX_CLASHES_1: Maximum allowable clashes per atom in the folded structure
# - global_score_weights: Weights applied to metrics to calculate the global ranking score
# - top_n_score_ESM: Number of top designs to keep based on the calculated ESM global score
# - top_n_gnina_ESM: Number of top designs to keep based on GNINA docking scores
clash_distance = config.get('clash_distance', 2.0)
bond_distance = config.get('bond_distance', 1.2)
MIN_PLDDT_1 = config.get('MIN_PLDDT_1', 0.8)
MIN_RMSD_1 = config.get('MIN_RMSD_1', 0.5)
MAX_RMSD_1 = config.get('MAX_RMSD_1', 10)
MAX_CLASHES_1 = config.get('MAX_CLASHES_1', 1.01)
global_score_weights = config.get('global_score_weights', {"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2})
top_n_score_ESM = config.get('top_n_score_ESM', 10)
top_n_gnina_ESM = config.get('top_n_gnina_ESM', 9)

### GNINA params
# - gnina_path: Command or path to execute gnina
# - gnina_cnn: CNN scoring model to use in gnina
# - gnina_exhaustiveness: Exhaustiveness of the global search
# - gnina_autobox_add: Padding added around the ligand to define the autobox for docking
gnina_path = config.get('gnina_path', "gnina")
gnina_cnn = config.get('gnina_cnn', "crossdock_default2018")
gnina_exhaustiveness = config.get('gnina_exhaustiveness', 8)
gnina_autobox_add = config.get('gnina_autobox_add', 4)

### Second prediction round
# - model_flag: Prediction model to use ("CHAI", "AF3", or "BOLTZ_ONLY")
# - msa_flag: If true, MSAs will be used during structure prediction
model_flag = config.get('model_flag', "CHAI")
msa_flag = config.get('msa_flag', False)

### Boltz params
# - max_dist: Maximum distance threshold for Boltz prediction
# - use_msa_boltz: If true, uses MSAs during Boltz prediction
# - use_forces: If true, physical forces are applied in Boltz
# - no_kernels: Flag to disable custom kernels in Boltz
# - path_to_boltz_env: Path to the environment needed to use Boltz
# - devices: Number of GPU devices to use for Boltz
# - recycling_steps: Number of recycling iterations for structure prediction
# - sampling_steps: Number of standard sampling steps
# - diffusion_samples: Number of diffusion samples to generate
# - output_format: Output format for the final structures ('pdb', 'cif', etc.)
# - sampling_steps_affinity: Number of sampling steps used specifically for affinity calculation
# - binding_pocket: List of residues defining the binding pocket
max_dist = config.get('max_dist', 5.0)
use_msa_boltz = config.get('use_msa_boltz', True)
use_forces = config.get('use_forces', True)
no_kernels = config.get('no_kernels', True)
path_to_boltz_env = config.get('path_to_boltz_env')
devices = config.get('devices', 1)
recycling_steps = config.get('recycling_steps', 3)
sampling_steps = config.get('sampling_steps', 100)
diffusion_samples = config.get('diffusion_samples', 1)
output_format = config.get('output_format', 'pdb')
sampling_steps_affinity = config.get('sampling_steps_affinity', 100)
binding_pocket = config.get('binding_pocket', None)

### Final filtering
# - top_n_score_final: Number of designs to be selected after the final round of filtering based on structural scores
# - top_n_gnina_final: Number of designs to be selected after the final round of filtering based on GNINA scores
top_n_score_final = config.get('top_n_score_final', 5)
top_n_gnina_final = config.get('top_n_gnina_final', 2)


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
if  RF_folder: # If RF_folder is provided, it means we are skipping RF diffusion and using existing designs, so we can skip ligand sampling if the sampled ligand pdbs already exist
    lowest_energy_conformer = func.generate_conformer(smiles=ligand_smiles, n=num_conformers, output_dir=f"{args.output}/ligand_sampling")
else:
    lowest_energy_conformer =run_ligand_sampling_pipeline(
            ligand_smiles=ligand_smiles,
            structure_path=structure_path,
            output_path=f"{args.output}/ligand_sampling",
            num_conformers=num_conformers,
            num_positions=num_positions,
            ligand_name=ligand_name, # Keep default or add CLI arg for this,
            rotation=rotation,
            translation=translation,
            conformer_rmsd_cutoff=conformer_rmsd_cutoff
            )

### Generate redesign information

pdb_info, printable_pdb_info = func.extract_pdb_info(pdb_path=structure_path,first_shell_distance=first_shell_distance,
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

if not ligand_MPNN_only or not RF_folder:
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
elif RF_folder:
    print(f"Using RF diffusion designs from {RF_folder} as input for Ligand MPNN...")
    mpnn_pdb_folder = RF_folder
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

func.run_ESMfold(f'{args.output}/MPNN_df.csv',f'{args.output}/ESM_predictions/',path_to_ESM_env=path_to_ESM_env,path_to_ESM_script=path_to_ESM_script,path_to_ESM_image=path_to_ESM_image)

### First 3D filter (3D quality)

ESM_df = func.threed_params_1_df(folder=f'{args.output}/ESM_predictions/', original_path=structure_path, clash_distance=clash_distance, bond_distance=bond_distance,
                                 ligand_path=lowest_energy_conformer, # I select lowest energy conformer as an assumption
                                 output_folder=f"{args.output}",
                                 output_name='ESM_filtered_df.csv',
                                 gnina_path=gnina_path,
                                 cnn=gnina_cnn, 
                                 exhaustiveness=gnina_exhaustiveness, 
                                 autobox_add=gnina_autobox_add)
ESM_df = func.threed_filter_1_df(ESM_df, output_folder=f"{args.output}", weights=global_score_weights, min_plddt=MIN_PLDDT_1, min_rmsd=MIN_RMSD_1, max_rmsd=MAX_RMSD_1, max_clashes=MAX_CLASHES_1, top_n_score=top_n_score_ESM,top_n_gnina=top_n_gnina_ESM)
ESM_df.to_csv(f'{args.output}/ESM_filtered_df.csv',index=False)

### Chai/AF prediction 

boltz_df = func.second_prediction_round(
    df=ESM_df, 
    model_flag=model_flag, 
    MSA_flag=msa_flag, 
    ligand_smiles=ligand_smiles, 
    output_path=f"{args.output}/{model_flag}_predictions",
    args_output=args.output, 
    # Boltz specific parameters
    pocket_list=None, 
    max_dist=5.0,
    use_msa_boltz=True,
    use_forces=True, 
    no_kernels=True,
    path_to_boltz_env="/home/eduardo/mambaforge/envs/boltz",
    devices=1, 
    recycling_steps=3, 
    sampling_steps=100, 
    diffusion_samples=1, 
    output_format='pdb', 
    sampling_steps_affinity=100
)
#boltz_df.to_csv(f'{args.output}/boltz_df.csv')
#If it's the first time running Chai, it will take longer as it needs to download the model weights

### Second 3D filter (3D quality and pocket metrics)
# Regular prediction
second_prediction_df = func.threed_params_1_df(folder=f'{args.output}/{model_flag}_pdbs/{model_flag}_prediction_pdbs', original_path=structure_path, clash_distance=clash_distance, bond_distance=bond_distance,
                                 ligand_path=lowest_energy_conformer, # I select lowest energy conformer as an assumption
                                 output_folder=f"{args.output}",
                                 output_name=f'{model_flag}_predictions.csv',
                                 gnina_path=gnina_path,
                                 cnn=gnina_cnn, 
                                 exhaustiveness=gnina_exhaustiveness, 
                                 autobox_add=gnina_autobox_add)
second_prediction_df.to_csv(f'{args.output}/{model_flag}_predictions.csv')

# Cofold prediction
second_prediction_cofold_df = func.threed_params_1_df(folder=f'{args.output}/{model_flag}_pdbs/{model_flag}_cofold_pdbs', original_path=structure_path, clash_distance=clash_distance, bond_distance=bond_distance,
                                 ligand_path=lowest_energy_conformer, # I select lowest energy conformer as an assumption
                                 output_folder=f"{args.output}",
                                 output_name=f'{model_flag}_cofold_predictions.csv',
                                 gnina_path=gnina_path,
                                 cnn=gnina_cnn, 
                                 exhaustiveness=gnina_exhaustiveness, 
                                 autobox_add=gnina_autobox_add)

second_prediction_cofold_df = func.check_cofold_validity(
    df=second_prediction_cofold_df, 
    path_to_structures=f'{args.output}/{model_flag}_pdbs/{model_flag}_cofold_pdbs', 
    ligand_id='LIG2', 
    site_residues=binding_pocket, 
    center_cutoff=4.0,
    extension=".cif")

second_prediction_cofold_df.to_csv(f'{args.output}/{model_flag}_cofold_predictions.csv')
# Boltz prediction
boltz_prediction_df = func.threed_params_1_df(folder=f'{args.output}/BOLTZ_pdbs', original_path=structure_path, clash_distance=clash_distance, bond_distance=bond_distance,
                                 ligand_path=lowest_energy_conformer, # I select lowest energy conformer as an assumption
                                 output_folder=f"{args.output}",
                                 output_name=f'BOLTZ_predictions.csv',
                                 gnina_path=gnina_path,
                                 cnn=gnina_cnn, 
                                 exhaustiveness=gnina_exhaustiveness, 
                                 autobox_add=gnina_autobox_add)

boltz_prediction_df = func.check_cofold_validity(
    df=second_prediction_cofold_df, 
    path_to_structures=f'{args.output}/BOLTZ_pdbs', 
    ligand_id='LIG', 
    site_residues=binding_pocket, 
    center_cutoff=4.0,
    extension=".pdb")

boltz_prediction_df.to_csv(f'{args.output}/BOLTZ_predictions.csv')


### Final filtering and selection of best candidates

func.threed_filter_2_df([second_prediction_df, second_prediction_cofold_df, boltz_prediction_df], output_folder=f"{args.output}",
                         output_names=[f'final{model_flag}_predictions.csv', f'final_{model_flag}_cofold_predictions.csv', f'final_BOLTZ_predictions.csv'],
                         weights=global_score_weights, min_plddt=MIN_PLDDT_1, min_rmsd=MIN_RMSD_1, max_rmsd=MAX_RMSD_1, max_clashes=MAX_CLASHES_1, top_n_score=top_n_score_final,top_n_gnina=top_n_gnina_final)



