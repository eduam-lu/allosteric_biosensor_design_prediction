import pandas as pd
import os
import subprocess
from pathlib import Path

# -----------------------------------------------------------------------------
# 1. SETUP: Scan the input folder
# -----------------------------------------------------------------------------
# glob_wildcards dynamically finds all PDBs in your starting folder.
# If there are 50 PDBs, SAMPLES will contain 50 unique names.
RF_DIR = "RF_final_pdbs"
SAMPLES, = glob_wildcards(f"{RF_DIR}/{{sample}}.pdb")

rule all:
    input:
        "results/final_top_3_designs.csv"

# -----------------------------------------------------------------------------
# 2. PHASE 1: SCATTER LigandMPNN (Parallel on GPUs)
# -----------------------------------------------------------------------------
rule run_ligand_mpnn:
    input:
        pdb = f"{RF_DIR}/{{sample}}.pdb"
    output:
        fasta = "MPNN_designs/{sample}_seqs.fa"
    resources: 
        gpu = 1
    conda: 
        "ligandmpnn_env"
    shell:
        """
        python /home/eduardo/LigandMPNN/run.py \
            --pdb_path {input.pdb} \
            --out_folder MPNN_designs \
            --batch_size 10
        """

# -----------------------------------------------------------------------------
# 3. PHASE 2: GATHER & BULK PREDICT (ESMfold)
# -----------------------------------------------------------------------------
# We gather all MPNN FASTAs, apply 1D filter, and run the bulk ESM script.
rule run_bulk_esmfold:
    input:
        fastas = expand("MPNN_designs/{sample}_seqs.fa", sample=SAMPLES)
    output:
        mpnn_csv = "results/MPNN_df.csv",
        esm_dir = directory("ESM_predictions")
    conda: 
        "esm_fold_env"
    run:
        from functions_bsd import process_MPNN_file, filter_dataframe_1D
        
        # Process individual FASTAs into one dataframe
        df_list = [process_MPNN_file(Path(f), top_n=5) for f in input.fastas]
        df = pd.concat(df_list, ignore_index=True)
        
        # 1D Filtering (Polyalanine/Polyglutamate)
        df = filter_dataframe_1D(df, window_size=10, threshold=10, aa='A')
        df = filter_dataframe_1D(df, window_size=10, threshold=10, aa='E')
        df.to_csv(output.mpnn_csv, index=False)
        
        # Run bulk ESM prediction
        Path(output.esm_dir).mkdir(parents=True, exist_ok=True)
        cmd = f"python esm_high_throughput.py --input_csv {output.mpnn_csv} --output_folder {output.esm_dir}"
        subprocess.run(cmd, shell=True, check=True)

# -----------------------------------------------------------------------------
# 4. PHASE 3: THE CHECKPOINT (3D Filter & GNINA)
# -----------------------------------------------------------------------------
# This is a 'checkpoint' instead of a 'rule'. Snakemake will re-evaluate 
# the rest of the pipeline AFTER this step finishes.
checkpoint filter_esm_structures:
    input:
        mpnn_csv = "results/MPNN_df.csv",
        esm_dir = "ESM_predictions"
    output:
        filtered_csv = "results/ESM_filtered_df.csv"
    run:
        from functions_bsd import threed_params_1_df, threed_filter_1_df
        
        # Generate 3D metrics & GNINA scores
        df_metrics = threed_params_1_df(
            folder=input.esm_dir, 
            output_folder="results",
            output_name="temp_metrics.csv",
            original_path="path/to/original.pdb" 
        )
        
        # Filter down to the survivors
        weights = {"pLDDT_mean": 0.4, "TMscore": 0.4, "clashes_per_atom": -0.2}
        df_filtered = threed_filter_1_df(df_metrics, output_folder="results", weights=weights)
        df_filtered.to_csv(output.filtered_csv, index=False)

# -----------------------------------------------------------------------------
# 5. PHASE 4: SCATTER ALPHAFOLD3 ON SURVIVORS
# -----------------------------------------------------------------------------
rule generate_af3_json:
    input: 
        csv = "results/ESM_filtered_df.csv"
    output: 
        json = "AF3_inputs/{id}.json"
    run:
        from functions_bsd import AF_json_generator_wo_MSA
        df = pd.read_csv(input.csv)
        row = df[df['file_ID'] == wildcards.id].iloc[0]
        AF_json_generator_wo_MSA(row, "AF3_inputs")
        # Rename the generic temp.json to the specific {id}.json
        Path("AF3_inputs/temp.json").rename(output.json)

rule run_alphafold3:
    input: 
        json = "AF3_inputs/{id}.json"
    output: 
        cif = "AF3_predictions/{id}/model.cif"
    resources: 
        gpu = 1
    conda: 
        "alphafold3_env"
    shell:
        """
        python /mnt/data/alphafold3/run_alphafold.py --json_path={input.json} --output_dir=AF3_predictions/{wildcards.id} ...
        """

# -----------------------------------------------------------------------------
# 6. PHASE 5: GATHER FINAL SURVIVORS
# -----------------------------------------------------------------------------
# This function dynamically looks at the checkpoint to see what IDs survived.
def get_surviving_af3_targets(wildcards):
    # This forces Snakemake to wait until 'filter_esm_structures' is done
    checkpoint_out = checkpoints.filter_esm_structures.get(**wildcards).output.filtered_csv
    
    # Read the surviving IDs
    df = pd.read_csv(checkpoint_out)
    surviving_ids = df['file_ID'].tolist()
    
    # Request an AF3 CIF file for every ID that survived
    return expand("AF3_predictions/{id}/model.cif", id=surviving_ids)

# The final rule uses the dynamic function as its input
rule final_selection:
    input:
        cifs = get_surviving_af3_targets
    output:
        "results/final_top_3_designs.csv"
    run:
        from functions_bsd import threed_params_1_df, threed_filter_2_df
        
        # Example logic: run metrics on the AF3 folder and select top 3
        # (You can use your actual threed_filter_2_df logic here)
        final_df = threed_params_1_df(folder="AF3_predictions", ...)
        final_top3 = final_df.sort_values(by="pLDDT_mean", ascending=False).head(3)
        final_top3.to_csv(output[0], index=False)