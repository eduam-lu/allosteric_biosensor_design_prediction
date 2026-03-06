# Strategy 1: Allosteric binding pocket redesign
## Strategy description

## Pipeline description

## Usage

### Config file

```yaml
### Ligand sampling params
# - ligand_smiles: SMILES string of the ligand to be sampled
# - ligand_name: Name of the ligand (used in PDB files)
# - num_conformers: Number of conformers to generate for the ligand
# - num_positions: Number of different positions to sample for the ligand around the protein
# - conformer_rmsd_cutoff: RMSD cutoff for filtering conformers. Anything below this value will be considered redundant
# - rotation: Rotation angle in degrees for exploring ligand placement
# - translation: Translation distance in angstroms for exploring ligand placement
ligand_smiles: "C[C@@]1([C@H]2C[C@H]3[C@@H](C(=O)C(=C([C@]3(C(=O)C2=C(C4=C1C=CC=C4O)O)O)O)C(=O)N)N(C)C)O"
ligand_name: "Tc"
num_conformers: 2
num_positions: 1
conformer_rmsd_cutoff: 0.75
rotation: 15
translation: 0.3

### Active site definition params
# - first_shell_distance: Distance in angstroms to define the first shell of residues around the ligand
# - second_shell_distance: Distance in angstroms to define the second shell of residues around the first shell
# - include_second_shell: If true, residues in the second shell will also be included in the redesign
# - user_defined_active_site: If provided, this list of residues will be used as the active site instead of calculating one
# - user_defined_residues: If provided, the program ensures this list of residues are included in the active site
first_shell_distance: 5.0
second_shell_distance: 5.0
include_second_shell: false
user_defined_active_site: null
user_defined_residues: null

### Redesign conditions params
# - ligand_MPNN_only: If true, only the ligandMPNN step will be performed, skipping RF diffusion
# - RF_folder: If provided, this folder with RF diffusion designs will be used as input for LigandMPNN instead of running RF diffusion
# - RF_model: RF diffusion model to use ("all_atom", "RF1", "RF3")
# - consertative_redesign: If true, residues within a secondary structure element will not be redesigned. Except if they are the first or final residue of the element (NOT IMPLEMENTED YET)
# - segment_extension: Number of residues to extend the redesign segment on each side
# - n_termini_extension: Number of residues to extend the redesign segment at the N-terminus
# - c_termini_extension: Number of residues to extend the redesign segment at the C-terminus. If null, extends to the end of the chain
# - user_defined_contig_map: If provided, this contig map will be used for RF diffusion redesign instead of calculating one
ligand_MPNN_only: false
RF_folder: null
RF_model: "all_atom"
consertative_redesign: false
segment_extension: 1
n_termini_extension: 1
c_termini_extension: null
user_defined_contig_map: null

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
path_to_RFAA_apptainer: "/home/eduardo/allostery/rf_diffusion_all_atom/rf_se3_diffusion.sif"
path_to_RFAA_script: "/home/eduardo/allostery/rf_diffusion_all_atom/run_inference.py"
path_to_RFAA_weights: "/home/eduardo/allostery/rf_diffusion_all_atom/RFDiffusionAA_paper_weights.pt"
path_to_RF1_script: "/home/eduardo/programs/RFdiffusion/scripts/run_inference.py"
path_to_RF1_env: "/home/sofia/miniforge3/envs/SE3-nvidia"
num_designs: 1
T: 50
RFAA_ligand_name: "UNL"
design_startnum: 1
deterministic: true

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
path_to_ligand_MPNN_script: "/home/eduardo/LigandMPNN/run.py"
path_to_ligand_MPNN_env: "/home/eduardo/miniforge3/envs/ligandmpnn_env"
mpnn_model_type: "ligand_mpnn"
path_to_mpnn_model: "/home/eduardo/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"
MPNN_num_designs: 10
n_batches: 1
mpnn_temperature: 0.05
bias_aa_global: null
omit_aa_global: null
side_chain_context: 0
first_shell_only: false
user_defined_mpnn_redesign: null
top_n_mpnn_candidates: 5

### ESM
# - path_to_ESM_env: Path to the ESM fold environment
# - path_to_ESM_script: Path to the high throughput ESM folding python script
# - path_to_ESM_image: Path to the ESM fold apptainer/singularity image
path_to_ESM_env: "/home/eduardo/miniforge3/envs/esm_fold_env"
path_to_ESM_script: "/home/eduardo/allostery/strategy_1/esm_high_throughput.py"
path_to_ESM_image: "/home/eduardo/conda_backups/esmfold_plus.sif"

### 1D filtering
# - filter_1d_window_size: Window size for detecting single-amino acid repeats (like Poly-A or Poly-E)
# - filter_1d_treshold: Threshold amount of repeating amino acids to trigger filtering a sequence
filter_1d_window_size: 10
filter_1d_treshold: 10

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
clash_distance: 2.0
bond_distance: 1.2
MIN_PLDDT_1: 0.8
MIN_RMSD_1: 0.5
MAX_RMSD_1: 10
MAX_CLASHES_1: 1.01
global_score_weights:
  pLDDT_mean: 0.4
  TMscore: 0.4
  clashes_per_atom: -0.2
top_n_score_ESM: 10
top_n_gnina_ESM: 9

### GNINA params
# - gnina_path: Command or path to execute gnina
# - gnina_cnn: CNN scoring model to use in gnina
# - gnina_exhaustiveness: Exhaustiveness of the global search
# - gnina_autobox_add: Padding added around the ligand to define the autobox for docking
gnina_path: "gnina"
gnina_cnn: "crossdock_default2018"
gnina_exhaustiveness: 8
gnina_autobox_add: 4

### Second prediction round
# - model_flag: Prediction model to use ("CHAI", "AF3", or "BOLTZ_ONLY")
# - msa_flag: If true, MSAs will be used during structure prediction
model_flag: "CHAI"
msa_flag: false

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
max_dist: 5.0
use_msa_boltz: true
use_forces: true
no_kernels: true
path_to_boltz_env: "/home/eduardo/mambaforge/envs/boltz"
devices: 1
recycling_steps: 3
sampling_steps: 100
diffusion_samples: 1
output_format: "pdb"
sampling_steps_affinity: 100
binding_pocket: null

### Final filtering
# - top_n_score_final: Number of designs to be selected after the final round of filtering based on structural scores
# - top_n_gnina_final: Number of designs to be selected after the final round of filtering based on GNINA scores
top_n_score_final: 5
top_n_gnina_final: 2```
