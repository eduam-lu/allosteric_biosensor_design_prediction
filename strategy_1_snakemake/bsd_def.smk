# -----------------------------------------------------------------------------
# IMPORTS
# -----------------------------------------------------------------------------
import os
import pandas as pd
from pathlib import Path

# Import your custom functions directly
import functions_bsd as func
# -----------------------------------------------------------------------------
# PARAMETERS
# -----------------------------------------------------------------------------

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
MPNN_GPUS = config.get("mpnn_gpus", 4)
SPLIT_IDS = [str(i) for i in range(MPNN_GPUS)]

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

# -----------------------------------------------------------------------------
# 1. SETUP: Scan the input RFdiffusion folder
# -----------------------------------------------------------------------------
# The inputs from the RFdiffusion module are a final pdb folder, containing one pdb per design. We use glob_wildcards to dynamically find all PDBs in your starting folder.
# Also, the module outputs the necessary MPNN json files, minimally pdb_paths_multijson and redesigned_residues_multijson
# Run this using: snakemake -s bsd_def.smk --configfile basic_config_bsd.yml --cores 8 --resources gpu=1

RF_PDB_DIR = config.get("rf_pdb_dir", "allostery/full_run/RF_final_pdbs")
OUTPUT_DIR = config.get("output_dir", "results")

# Dynamically find all PDBs in your starting folder
RF_DESIGNS, = glob_wildcards(os.path.join(RF_PDB_DIR, "{rf_id}.pdb"))

# Target rule to drive the DAG
rule all:
    input:
        f"{OUTPUT_DIR}/final_scores_topn_filtered.csv"
# -----------------------------------------------------------------------------
# 2. PHASE 1: SCATTER LigandMPNN (Parallel on GPUs)
# -----------------------------------------------------------------------------
#Here a LigandMPNN rule will take each pdb, run the MPNN, and output a fasta file for each pdb. This is a scatter step that can be fully parallelized on GPUs.
# Then in the processing rule, the MPNN outputs will be read, processed into a dataframe, and filtered by MPNN score and by the 1D filter. The resulting dataframe will be saved as a csv for the next step.

# I think the best way to parallelise is a hybrid strategy
# --- A. Split the JSONs ---
# This rule takes the JSONs from RFdiffusion and divides them equally 
# into N separate chunks (where N = MPNN_GPUS).
rule split_mpnn_jsons:
    input:
        # Assuming these are generated by your RFdiffusion module beforehand
        paths_json = config.get("rf_paths_json", f"{OUTPUT_DIR}/RF_outputs/pdb_paths_multi.json"),
        res_json = config.get("rf_res_json", f"{OUTPUT_DIR}/RF_outputs/redesigned_residues_multi.json")
    output:
        # We output a list of N json files
        paths_splits = expand(f"{OUTPUT_DIR}/MPNN_jsons/split_{{split_id}}/pdb_paths_multi.json", split_id=SPLIT_IDS),
        res_splits = expand(f"{OUTPUT_DIR}/MPNN_jsons/split_{{split_id}}/redesigned_residues_multi.json", split_id=SPLIT_IDS)
    run:
        import json
        import math
        
        # Load the original JSONs
        with open(input.paths_json, 'r') as f: paths = json.load(f)
        with open(input.res_json, 'r') as f: res = json.load(f)
            
        keys = list(paths.keys())
        n_splits = len(output.paths_splits)
        
        # Calculate how many PDBs go into each chunk
        chunk_size = math.ceil(len(keys) / n_splits) if keys else 1
        
        for i in range(n_splits):
            start_idx = i * chunk_size
            end_idx = start_idx + chunk_size
            chunk_keys = keys[start_idx:end_idx]
            
            # Slice the dictionaries
            chunk_paths = {k: paths[k] for k in chunk_keys}
            chunk_res = {k: res.get(k, "") for k in chunk_keys}
            
            # Create directories and write the split JSONs
            os.makedirs(os.path.dirname(output.paths_splits[i]), exist_ok=True)
            with open(output.paths_splits[i], 'w') as f: json.dump(chunk_paths, f)
            with open(output.res_splits[i], 'w') as f: json.dump(chunk_res, f)

# --- B. Run LigandMPNN on each split ---
# Snakemake will automatically launch this rule N times concurrently, 
# assigning 1 GPU to each chunk.
rule run_ligand_mpnn_hybrid:
    input:
        paths_json = f"{OUTPUT_DIR}/MPNN_jsons/split_{{split_id}}/pdb_paths_multi.json",
        res_json = f"{OUTPUT_DIR}/MPNN_jsons/split_{{split_id}}/redesigned_residues_multi.json"
    output:
        seqs_dir = directory(f"{OUTPUT_DIR}/mpnn_outputs/split_{{split_id}}/seqs")
    params:
        out_folder = f"{OUTPUT_DIR}/mpnn_outputs/split_{{split_id}}"
    resources:
        gpu = 1
    run:
        from snakemake.shell import shell
        import json
        
        # Check if chunk is empty
        with open(input.paths_json, 'r') as f:
            has_keys = len(json.load(f)) > 0
            
        if has_keys:
            # Base command
            cmd_parts = [
                f"conda run -p {config['path_to_ligand_MPNN_env']} python {config['path_to_ligand_MPNN_script']}",
                f"--pdb_path_multi {input.paths_json}",
                f"--out_folder {params.out_folder}",
                f"--save_stats 1",
                f"--batch_size {config.get('MPNN_num_designs', 10)}",
                f"--number_of_batches {config.get('n_batches', 1)}",
                f"--model_type {config.get('mpnn_model_type', 'ligand_mpnn')}",
                f"--checkpoint_ligand_mpnn {config['path_to_mpnn_model']}",
                f"--temperature {config.get('mpnn_temperature', 0.05)}",
                f"--ligand_mpnn_use_side_chain_context {config.get('side_chain_context', 0)}",
                f"--redesigned_residues_multi \"{input.res_json}\""
            ]
            
            # Connect the optional bias/omit parameters!
            if config.get("omit_aa_global"):
                cmd_parts.append(f"--omit_AA {config['omit_aa_global']}")
            if config.get("bias_aa_global"):
                cmd_parts.append(f"--bias_AA {config['bias_aa_global']}")
                
            shell(" ".join(cmd_parts))
        else:
            print(f"Chunk {wildcards.split_id} is empty. Skipping LigandMPNN execution.")
            shell(f"mkdir -p {output.seqs_dir}")

# --- C. Gather, Process, and Filter ---
# This rule waits for ALL splits to finish, gathers all the seqs/ directories,
# merges them into a single dataframe, and applies your 1D filters.
rule process_mpnn_output_and_filter:
    input:
        seqs_dirs = expand(f"{OUTPUT_DIR}/mpnn_outputs/split_{{split_id}}/seqs", split_id=SPLIT_IDS)
    output:
        mpnn_filtered_csv = f"{OUTPUT_DIR}/MPNN_filtered.csv"
    run:
        all_dfs = []
        
        # Iterate over all the seq directories generated by our hybrid approach
        for seq_dir in input.seqs_dirs:
            if os.path.exists(seq_dir) and any(os.scandir(seq_dir)):
                df = func.process_MPNN_folder(
                    folder=seq_dir,
                    top_n=config.get("top_n_mpnn_candidates", 5)
                )
                all_dfs.append(df)
        
        if not all_dfs:
            raise ValueError("No MPNN sequences were generated across any splits!")
            
        # Merge all splits back into one main dataframe
        mpnn_df = pd.concat(all_dfs, ignore_index=True)
        
        # Apply your 1D filters
        win_size = config.get("filter_1d_window_size", 10)
        thresh = config.get("filter_1d_treshold", 10)
        
        mpnn_df = func.filter_dataframe_1D(mpnn_df, window_size=win_size, threshold=thresh, seq_col='seq', aa='A')
        mpnn_df = func.filter_dataframe_1D(mpnn_df, window_size=win_size, threshold=thresh, seq_col='seq', aa='E')
        
        mpnn_df.to_csv(output.mpnn_filtered_csv, index=False)

# -----------------------------------------------------------------------------
# 3. PHASE 2: SCATTER ESMfold (Parallel on GPUs)
# -----------------------------------------------------------------------------
# Here, we will run ESMfold on all the MPNN designs that passed the 1D filter. This is another scatter step that can be fully parallelized on GPUs. The outputs will be the predicted structures for each design.
# So we need to adapt the ESM_highthroughput script to be scatter rather than bulk. Each ESMfold run will take one sequence, predict its structure, and save the output pdb.
# I am worried however, because the model takes time to load, so i wish there was a way to just load it once and then run all predictions sequentially. 
# If not, we will have to stick to bulk highthroughput.

# --- A. Split the MPNN output CSV ---
rule split_esm_csv:
    input:
        mpnn_filtered_csv = rules.process_mpnn_output_and_filter.output.mpnn_filtered_csv
    output:
        split_csvs = expand(f"{OUTPUT_DIR}/ESM_predictions/split_{{split_id}}/input.csv", split_id=ESM_SPLIT_IDS)
    run:
        import numpy as np
        
        df = pd.read_csv(input.mpnn_filtered_csv)
        n_splits = len(output.split_csvs)
        
        # np.array_split handles unequal divisions gracefully 
        # (e.g., splitting 10 rows into 4 chunks = 3, 3, 2, 2)
        df_splits = np.array_split(df, n_splits)
        
        for i, split_df in enumerate(df_splits):
            os.makedirs(os.path.dirname(output.split_csvs[i]), exist_ok=True)
            split_df.to_csv(output.split_csvs[i], index=False)

# --- B. Run ESMfold on each split ---
rule run_esmfold_hybrid:
    input:
        split_csv = f"{OUTPUT_DIR}/ESM_predictions/split_{{split_id}}/input.csv"
    output:
        esm_out_dir = directory(f"{OUTPUT_DIR}/ESM_predictions/split_{{split_id}}/pdbs")
    resources:
        gpu = 1
    conda:
        config[path_to_ESM_env]
    shell:
        """
        LINES=$(wc -l < {input.split_csv})
        if [ "$LINES" -gt "1" ]; then
            # Connected the python script path properly here:
             python {config[path_to_ESM_script]} \
                --input_csv {input.split_csv} \
                --output_folder {output.esm_out_dir}
        else
            mkdir -p {output.esm_out_dir}
        fi
        """

# --- C. Gather the PDBs back into one folder ---
# This ensures that the downstream checkpoint doesn't need to be rewritten 
# to loop through multiple directories.
rule gather_esm_pdbs:
    input:
        esm_dirs = expand(f"{OUTPUT_DIR}/ESM_predictions/split_{{split_id}}/pdbs", split_id=ESM_SPLIT_IDS)
    output:
        gathered_dir = directory(f"{OUTPUT_DIR}/ESM_predictions/all_pdbs")
    run:
        import shutil
        os.makedirs(output.gathered_dir, exist_ok=True)
        
        for d in input.esm_dirs:
            if os.path.exists(d):
                for f in os.listdir(d):
                    if f.endswith(".pdb"):
                        shutil.copy2(os.path.join(d, f), os.path.join(output.gathered_dir, f))


# -----------------------------------------------------------------------------
# 4. PHASE 3: FIRST CHECKPOINT (3D Params, Filter & GNINA)
# -----------------------------------------------------------------------------
# This is a 'checkpoint' instead of a 'rule'. Snakemake will re-evaluate 
# the rest of the pipeline AFTER this step finishes.

# First, we will calculate the 3D parameters for all ESMfold predictions. Then, we will apply a 3D filter to select survivors. 
#Finally, we will run GNINA on the survivors to get binding scores and filter by that as well. 
# The final output of this checkpoint will be a csv with the survivors and their scores, which will be the input for the next phase.

checkpoint first_3d_filter:
    input:
        esm_dir = rules.gather_esm_pdbs.output.gathered_dir
    output:
        survivors_csv = f"{OUTPUT_DIR}/checkpoints/ESM_survivors.csv"
    run:
        esm_df = func.threed_params_1_df(
            folder=input.esm_dir, 
            output_folder=f"{OUTPUT_DIR}/checkpoints",
            output_name="ESM_metrics.csv", 
            original_path=config["structure_path"], 
            # Reconnected all the 3D and GNINA parameters!
            clash_distance=config.get("clash_distance", 2.0), 
            bond_distance=config.get("bond_distance", 1.2),
            ligand_path=config["ligand_smiles"], 
            gnina_path=config.get("gnina_path", "gnina"),
            cnn=config.get("gnina_cnn", "crossdock_default2018"),
            exhaustiveness=config.get("gnina_exhaustiveness", 8),
            autobox_add=config.get("gnina_autobox_add", 4)
        )
        
        filtered_df = func.threed_filter_1_df(
            df=esm_df, 
            output_folder=f"{OUTPUT_DIR}/checkpoints", 
            weights=config.get("global_score_weights", {"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2}), 
            min_plddt=config.get("MIN_PLDDT_1", 0.8), 
            min_rmsd=config.get("MIN_RMSD_1", 0.5), 
            max_rmsd=config.get("MAX_RMSD_1", 10), 
            max_clashes=config.get("MAX_CLASHES_1", 1.01), 
            top_n_score=config.get("top_n_score_ESM", 10),
            top_n_gnina=config.get("top_n_gnina_ESM", 9)
        )
        filtered_df.to_csv(output.survivors_csv, index=False)


# Helper function to dynamically read the checkpoint output
def get_survivor_ids(wildcards):
    ckpt_out = checkpoints.first_3d_filter.get(**wildcards).output.survivors_csv
    df = pd.read_csv(ckpt_out)
    # Assumes 'file_ID' holds the design identifier (e.g., 'design_1')
    return df['file_ID'].tolist()

# -----------------------------------------------------------------------------
# 5. PHASE 4: SCATTER Boltz ON SURVIVORS (Parallel on GPUs)
# -----------------------------------------------------------------------------

# I will run each of the surviving sequences through Boltz, in an scatter manner. Keep in mind that each sequence requires one yaml file with the Boltz parameters, which will be generated from the csv of survivors. 
# The outputs will be the Boltz-designed structures for each survivo.
rule prep_boltz_yaml:
    input:
        survivors_csv = f"{OUTPUT_DIR}/checkpoints/ESM_survivors.csv"
    output:
        yaml = f"{OUTPUT_DIR}/boltz_yamls/{{survivor_id}}.yaml"
    run:
        df = pd.read_csv(input.survivors_csv)
        row = df[df['file_ID'] == wildcards.survivor_id].iloc[0]
        
        # Use your custom generator
        func.boltz_yaml_generator(
            row=row, 
            yaml_path=f"{OUTPUT_DIR}/boltz_yamls", 
            ligand_smiles=config["ligand_smiles"], 
            pocket_list=config.get("binding_pocket", []), 
            max_dist=config["max_dist"]
        )

#rule run_boltz:
rule run_boltz:
    input:
        yaml = rules.prep_boltz_yaml.output.yaml
    output:
        # A hidden flag file to tell Snakemake this specific survivor is done
        done_flag = f"{OUTPUT_DIR}/BOLTZ_raw/{{survivor_id}}/.done"
    resources:
        gpu = config.get("devices", 1)
    params:
        out_dir = f"{OUTPUT_DIR}/BOLTZ_raw",
        # Pull all Boltz parameters from your basic_config_bsd.yml
        env = config.get("path_to_boltz_env"),
        devices = config.get("devices", 1),
        recycling = config.get("recycling_steps", 3),
        sampling = config.get("sampling_steps", 100),
        diffusion = config.get("diffusion_samples", 1),
        out_format = config.get("output_format", "pdb"),
        affinity_sampling = config.get("sampling_steps_affinity", 100),
        use_msa = config.get("use_msa_boltz", True),
        use_forces = config.get("use_forces", True),
        no_kernels = config.get("no_kernels", True)
    run:
        from snakemake.shell import shell
        import os

        # 1. Base command with mandatory flags
        # Note: Using '-p {params.env}' to safely use the explicit path from your config
        command_parts = [
            f"conda run -p {params.env} boltz predict",
            f"{input.yaml}",
            f"--out_dir={params.out_dir}",
            f"--devices={params.devices}",
            f"--recycling_steps={params.recycling}",
            f"--sampling_steps={params.sampling}",
            f"--diffusion_samples={params.diffusion}",
            f"--output_format={params.out_format}",
            f"--sampling_steps_affinity={params.affinity_sampling}"
        ]

        # 2. Conditional flags
        if params.use_msa:
            command_parts.append("--use_msa_server")
        
        if params.use_forces:
            command_parts.append("--use_potentials")
            
        if params.no_kernels:
            command_parts.append("--no_kernels")

        # 3. Join all parts into a single string
        boltz_command = " ".join(command_parts)
        
        # 4. Run the process (Snakemake's shell automatically handles errors/logs)
        shell(boltz_command)

        # 5. Create the flag file so the process_boltz_folder rule knows it finished
        os.makedirs(f"{params.out_dir}/{wildcards.survivor_id}", exist_ok=True)
        with open(output.done_flag, "w") as f:
            f.write("done")

        # 6. Remove yaml file (as per the comment in your function!)
        os.remove(input.yaml)

# rule process_boltz_folder:
rule process_boltz_folder:
    input:
        # We use 'expand' to wait for ALL scattered Boltz jobs to finish.
        # Assuming run_boltz outputs a dummy file or specific nested file to mark completion:
        boltz_dones = expand(f"{OUTPUT_DIR}/BOLTZ_raw/{{survivor_id}}/.done", survivor_id=get_survivor_ids)
    output:
        # The flattened directory where the function will copy the PDBs
        pdbs_dir = directory(f"{OUTPUT_DIR}/BOLTZ_pdbs"),
        # The dataframe containing confidence, iptm, affinity, etc.
        metrics_csv = f"{OUTPUT_DIR}/BOLTZ_native_metrics.csv"
    params:
        # The base directory where Boltz generated its nested outputs
        raw_dir = f"{OUTPUT_DIR}/BOLTZ_raw"
    run:
        # 1. Call your custom function to parse the raw Boltz outputs
        boltz_df = func.process_Boltz_folder(
            input_folder=params.raw_dir,
            output_pdbs_folder=output.pdbs_dir
        )
        
        # 2. Save the resulting dataframe to a CSV
        if not boltz_df.empty:
            boltz_df.to_csv(output.metrics_csv, index=False)
        else:
            print("Warning: process_Boltz_folder returned an empty dataframe!")
            # Create an empty file to satisfy Snakemake's output requirement
            open(output.metrics_csv, 'a').close()

# -----------------------------------------------------------------------------
# 6. PHASE 5: SCATTER SECOND PREDICTION (Parallel on GPUs)
# -----------------------------------------------------------------------------
# We will also run each of the surviving sequences through either Chai or AF3, chosen by the user. This is another scatter step that can be parallelized on GPUs.
# For each surviving sequence, a regular and a cofold prediction will be run. Keep in mind that also the user will configure wether MSAs and or templates are used 
#The outputs will be the predicted structures and cofold structures for each survivor after Boltz design.

rule run_second_prediction:
    input:
        survivors_csv = f"{OUTPUT_DIR}/checkpoints/ESM_survivors.csv"
    output:
        # A hidden flag file to tell Snakemake this specific survivor is done
        done_flag = f"{OUTPUT_DIR}/{config['model_flag']}_raw/{{survivor_id}}/.done"
    resources:
        gpu = 1
    run:
        df = pd.read_csv(input.survivors_csv)
        row = df[df['file_ID'] == wildcards.survivor_id].iloc[0]
        
        model = config["model_flag"]
        msa = config["msa_flag"]
        
        # Run the appropriate model into the base OUTPUT_DIR 
        # (Your functions automatically create 'CHAI_prediction', 'CHAI_cofold', etc. inside it)
        if model == "CHAI":
            if msa:
                func.run_chai_w_MSA(row, f"{OUTPUT_DIR}/{model}_prediction/{wildcards.survivor_id}")
                func.run_chai_cofold_w_MSA(row, f"{OUTPUT_DIR}/{model}_cofold/{wildcards.survivor_id}", config["ligand_smiles"])
            else:
                func.run_chai(row, f"{OUTPUT_DIR}/{model}_prediction/{wildcards.survivor_id}")
                func.run_chai_cofold(row, f"{OUTPUT_DIR}/{model}_cofold/{wildcards.survivor_id}", config["ligand_smiles"])
        elif model == "AF3":
            func.run_AlphaFold3(row, f"{OUTPUT_DIR}", msa_flag=msa, ligand_SMILES=config["ligand_smiles"])
            
        # Create the flag file so the processing rule knows this sequence finished
        import os
        os.makedirs(os.path.dirname(output.done_flag), exist_ok=True)
        with open(output.done_flag, "w") as f:
            f.write("done")

# rule process_second_prediction_folder:
rule process_second_prediction_folder:
    input:
        # Wait for all scattered Chai/AF3 jobs to finish by looking for the .done flags
        dones = expand(f"{OUTPUT_DIR}/{config['model_flag']}_raw/{{survivor_id}}/.done", survivor_id=get_survivor_ids)
    output:
        # The custom functions create these two specific subfolders inside the _pdbs directory
        pred_pdbs = directory(f"{OUTPUT_DIR}/{config['model_flag']}_pdbs/{config['model_flag']}_prediction_pdbs"),
        cofold_pdbs = directory(f"{OUTPUT_DIR}/{config['model_flag']}_pdbs/{config['model_flag']}_cofold_pdbs")
    params:
        # The base path where the raw folders (CHAI_prediction / CHAI_cofold) are located
        raw_base = f"{OUTPUT_DIR}",
        # The base path where the cleaned PDBs/CIFs will go
        out_base = f"{OUTPUT_DIR}/{config['model_flag']}_pdbs",
        model = config["model_flag"]
    run:
        if params.model == "AF3":
            func.process_AF3_folder(params.raw_base, params.out_base)
        elif params.model == "CHAI":
            func.process_chai_folder(params.raw_base, params.out_base)
        else:
            raise ValueError(f"Unknown model_flag: {params.model}")
# -----------------------------------------------------------------------------
# 7. PHASE 6: GATHER FINAL CANDIDATES (3D Params, Filter & GNINA)
# -----------------------------------------------------------------------------
# At this step we have three kinds of structures for each survivor: Boltz, regular prediction, and cofold prediction. 
# We will calculate 3D parameters for all of them, apply a final 3D filter, and run GNINA to get the final scores. 
# We will get a csv file for each, unfiltered, but already including the global scores
# Finally, for the final selection, the three global scores will be averaged to get a final score for each design, and the top n designs will be selected as the final output of the pipeline.

rule gather_final_candidates:
    input:
        boltz_pdbs = lambda w: expand(f"{OUTPUT_DIR}/BOLTZ_pdbs/{{survivor_id}}.pdb", survivor_id=get_survivor_ids(w)),
        reg_pdbs = rules.process_second_prediction_folder.output.pred_pdbs,
        cofold_pdbs = rules.process_second_prediction_folder.output.cofold_pdbs
    output:
        final_summary = f"{OUTPUT_DIR}/final_scores_topn_filtered.csv"
    run:
        model = config["model_flag"]
        
        # 1. Calculate Parameters for Boltz
        boltz_df = func.threed_params_1_df(
            folder=f"{OUTPUT_DIR}/BOLTZ_pdbs", 
            output_folder=OUTPUT_DIR,
            output_name="BOLTZ_predictions.csv",
            original_path=config["structure_path"],
            ligand_path=config["ligand_smiles"] # Assuming you parse this dynamically later
        )
        # Check validity for Boltz
        boltz_df = func.check_cofold_validity(boltz_df, f"{OUTPUT_DIR}/BOLTZ_pdbs", "LIG", config.get("binding_pocket"), extension=".pdb")
        
        # 2. Calculate Parameters for Regular Prediction
        reg_df = func.threed_params_1_df(
            folder=f"{OUTPUT_DIR}/{model}_pdbs", # Depends on how run_chai structures output
            output_folder=OUTPUT_DIR,
            output_name=f"{model}_predictions.csv",
            original_path=config["structure_path"]
        )

        # 3. Calculate Parameters for Cofold
        cofold_df = func.threed_params_1_df(
            folder=f"{OUTPUT_DIR}/{model}_pdbs", 
            output_folder=OUTPUT_DIR,
            output_name=f"{model}_cofold_predictions.csv",
            original_path=config["structure_path"],
            ligand_path=config["ligand_smiles"]
        )
        cofold_df = func.check_cofold_validity(cofold_df, f"{OUTPUT_DIR}/{model}_pdbs", "LIG2", config.get("binding_pocket"), extension=".cif")

        # 4. Final selection merging all three using your multi-filter
        weights = {"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2}
        
        func.threed_filter_2_df(
            df_list=[reg_df, cofold_df, boltz_df], 
            output_folder=OUTPUT_DIR,
            output_names=[f'final_{model}_predictions.csv', f'final_{model}_cofold.csv', 'final_BOLTZ.csv'],
            weights=config.get("global_score_weights", {"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2}), 
            min_plddt=config.get("MIN_PLDDT_1", 0.8), 
            min_rmsd=config.get("MIN_RMSD_1", 0.5), 
            max_rmsd=config.get("MAX_RMSD_1", 10), 
            max_clashes=config.get("MAX_CLASHES_1", 1.01), 
            top_n_score=config.get("top_n_score_final", 5),   
            top_n_gnina=config.get("top_n_gnina_final", 2)    
        )

