# ============================================================================
# RF DIFFUSION BINDING SITE DESIGNER - SNAKEMAKE PIPELINE
# ============================================================================
# Pipeline for binding site redesign using RF Diffusion and LigandMPNN
# This pipeline mimics the logic from binding_site_designer.py up to ligand MPNN
#
# Usage:
#   snakemake --snakefile rf_diffusion_bsd_snakemake.smk -c 1 --config config=config_rfdiffusion_bsd.yml
#
# ============================================================================
# IMPORTS
# ============================================================================

import os
import sys
import json
import yaml
import pandas as pd
from pathlib import Path
import shutil
import zipfile
import gzip
from Bio.PDB import PDBParser, PDBIO
from snakemake.shell import shell

# Import functions from functions_bsd.py
sys.path.insert(0, os.path.dirname(__file__))
import functions_snakemake_rfd as func

# ============================================================================
# CONFIGURATION LOADING
# ============================================================================

# Load configuration file
configfile: "config_rfdiffusion_bsd.yml"

# ============================================================================
# PARAMETERS - All config values loaded with defaults
# ============================================================================

### Basic parameters
structure_path = config.get('structure_path')
chain_id = config.get('chain_id', None)  # If None, will use first chain with residues
ligand_resname = config.get('ligand_resname', None)
ligand_name = config.get('ligand_name', "Tc")
ligand_smiles = config.get('ligand_smiles')
output_dir = config.get('output_dir', "outputs/rf_diffusion_bsd")

### Ligand sampling parameters
num_conformers = config.get('num_conformers', 2)
num_positions = config.get('num_positions', 1)
conformer_rmsd_cutoff = config.get('conformer_rmsd_cutoff', 0.75)
rotation = config.get('rotation', 15)
translation = config.get('translation', 0.3)

### Active site definition parameters
first_shell_distance = config.get('first_shell_distance', 5.0)
second_shell_distance = config.get('second_shell_distance', 5.0)
include_second_shell = config.get('include_second_shell', False)
user_defined_active_site = config.get('user_defined_active_site', None)
user_defined_residues = config.get('user_defined_residues', None)

### Redesign conditions parameters
ligand_MPNN_only = config.get('ligand_MPNN_only', False)
RF_folder = config.get('RF_folder', None)
RF_model = config.get('RF_model', "all_atom")
conservative_redesign = config.get('conservative_redesign', False)
segment_extension = config.get('segment_extension', 1)
n_termini_extension = config.get('n_termini_extension', 1)
c_termini_extension = config.get('c_termini_extension', None)
user_defined_contig_map = config.get('user_defined_contig_map', None)


### RF diffusion execution parameters
path_to_RFAA_apptainer = config.get('path_to_RFAA_apptainer')
path_to_RFAA_script = config.get('path_to_RFAA_script')
path_to_RFAA_weights = config.get('path_to_RFAA_weights')
path_to_RF1_script = config.get('path_to_RF1_script')
path_to_RF1_env = config.get('path_to_RF1_env')
num_designs = config.get('num_designs', 10)
T = config.get('T', 50)
RFAA_ligand_name = config.get('RFAA_ligand_name', "UNL")
design_startnum = config.get('design_startnum', 1)
deterministic = config.get('deterministic', True)

### RF3 specific parameters
path_to_RF3_env = config.get('path_to_RF3_env')
RF3_batch_size = config.get('RF3_batch_size', 1)
RF3_num_batches = config.get('RF3_num_batches', num_designs)
RF3_num_timesteps = config.get('RF3_num_timesteps', 200)
RF3_step_scale = config.get('RF3_step_scale', 1.5)
RF3_noise_scale = config.get('RF3_noise_scale', 1.003)
RF3_gamma_0 = config.get('RF3_gamma_0', 0.6)
RF3_gamma_min = config.get('RF3_gamma_min', 1.0)
RF3_cfg_scale = config.get('RF3_cfg_scale', 1.5)
RF3_use_classifier_free_guidance = config.get('RF3_use_classifier_free_guidance', False)
RF3_dump_trajectories = config.get('RF3_dump_trajectories', False)
RF3_prevalidate_inputs = config.get('RF3_prevalidate_inputs', False)
RF3_low_memory_mode = config.get('RF3_low_memory_mode', False)
RF3_checkpoint_path = config.get('RF3_checkpoint_path', 'rfd3')


# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def get_prepared_pdb():
    """Get the prepared PDB file after cleaning."""
    return f"{output_dir}/PDB_prep/cleaned_structure.pdb"

def get_pdb_info_file():
    """Get the PDB info JSON file."""
    return f"{output_dir}/PDB_info/pdb_info.json"

def get_lowest_energy_conformer():
    """Get the lowest energy ligand conformer."""
    ligand_sampling_dir = f"{output_dir}/ligand_sampling"
    return f"{ligand_sampling_dir}/lowest_energy_conformer.pdb"

def get_rf_model():
    """Get the RF model to use."""
    return RF_model

# ============================================================================
# 0. RULE ALL - TARGET RULE DRIVING THE DAG
# ============================================================================

rule all:
    """
    Target rule to generate RF Diffusion designs.
    This is the final output that drives the entire pipeline.
    """
    input:
        rf_designs_done = f"{output_dir}/RF_designs/designs_completed.txt",
        rf_final_cif_done = f"{output_dir}/RF_final_cif/unzip_completed.txt"
    message:
        "RF Diffusion binding site design pipeline completed successfully!"

# ============================================================================
# 1. CLEAN AND PREPARE PDB FILE
# ============================================================================

rule clean_pdb:
    """
    Remove water and other non-protein and non-ligand atoms.
    Keep only protein residues and specified ligands.
    Ensure only a single chain remains.
    Prepare the PDB file for input to RF Diffusion.
    """
    input:
        pdb = structure_path
    output:
        cleaned_pdb = f"{output_dir}/PDB_prep/cleaned_structure.pdb"
    params:
        output_dir = f"{output_dir}/PDB_prep",
        chain_id = chain_id,
        ligand_resname = ligand_resname
    run:
        
        Path(params.output_dir).mkdir(parents=True, exist_ok=True)
        
        # Parse the input PDB
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('protein', input.pdb)
        
        # Parse ligand residue names to keep (comma/space separated or single value)
        ligands_to_keep = set()
        if params.ligand_resname:
            # Handle both comma-separated and space-separated values
            ligand_list = str(params.ligand_resname).replace(',', ' ').split()
            ligands_to_keep = set([lig.strip().upper() for lig in ligand_list if lig.strip()])
        
        # Iterate through structure and handle chains
        model = structure[0]  # Get first model
        chains_to_remove = []
        selected_chain_id = params.chain_id
        
        # If chain_id not specified, select the first chain with protein residues
        if not selected_chain_id:
            for chain in model:
                has_protein = any(res.get_id()[0] == ' ' for res in chain)
                if has_protein:
                    selected_chain_id = chain.get_id()
                    break
        
        # Identify chains to remove (keep only the selected chain)
        for chain in model:
            if chain.get_id() != selected_chain_id:
                chains_to_remove.append(chain.get_id())
        
        # Remove unwanted chains
        for chain_id_to_remove in chains_to_remove:
            model.detach_child(chain_id_to_remove)
        
        print(f"✓ Keeping chain {selected_chain_id}")
        if chains_to_remove:
            print(f"  Removed chains: {', '.join(chains_to_remove)}")
        
        # Clean residues in the selected chain
        chain = model[selected_chain_id]
        residues_to_remove = []
        for residue in chain:
            hetflag = residue.get_id()[0]
            resname = residue.get_resname().strip().upper()
            
            # Mark for removal if heteroatom
            if hetflag != ' ':
                # Always remove water (HOH, WAT, etc.)
                if resname == 'HOH' or resname == 'WAT':
                    residues_to_remove.append(residue)
                # Keep only if it's in the ligands_to_keep list
                elif ligands_to_keep and resname not in ligands_to_keep:
                    residues_to_remove.append(residue)
                # If no ligands specified, keep all non-water heteroatoms
                elif not ligands_to_keep:
                    residues_to_remove.append(residue)
        
        # Remove marked residues
        for residue in residues_to_remove:
            chain.detach_child(residue.get_id())
        
        # Write cleaned structure
        io = PDBIO()
        io.set_structure(structure)
        io.save(output.cleaned_pdb)
        
        shell("echo 'PDB cleaning completed' > {params.output_dir}/clean_status.txt")

# ============================================================================
# 2. EXTRACT PDB INFORMATION
# ============================================================================

rule extract_pdb_info:
    """
    Extract PDB information including:
    - Chain ID
    - Sequence (SEQRES and coordinate-based)
    - Secondary structure (DSSP)
    - Active site definition
    
    Produces a JSON file with all metadata.
    """
    input:
        pdb = get_prepared_pdb()
    output:
        pdb_info = f"{output_dir}/PDB_info/pdb_info.json"
    params:
        output_dir = f"{output_dir}/PDB_info",
        first_shell_distance = first_shell_distance,
        second_shell_distance = second_shell_distance,
        user_defined_active_site = user_defined_active_site,
        user_defined_residues = user_defined_residues
    run:
        Path(params.output_dir).mkdir(parents=True, exist_ok=True)
        
        # Extract PDB information using functions_bsd
        pdb_info, printable_pdb_info = func.extract_pdb_info(
            pdb_path=input.pdb,
            first_shell_distance=params.first_shell_distance,
            second_shell_distance=params.second_shell_distance,
            user_defined_active_site=params.user_defined_active_site,
            user_defined_residues=params.user_defined_residues
        )
        
        # Save PDB info to JSON
        with open(output.pdb_info, 'w') as f:
            json.dump(printable_pdb_info, f, indent=2)
        
        print(f" PDB info extracted and saved to {output.pdb_info}")

# ============================================================================
# 3. LIGAND SAMPLING
# ============================================================================

rule ligand_sampling:
    """
    Sample ligand positions and conformers around the protein.
    Generates multiple conformers and positions within the binding pocket.
    """
    input:
        pdb = get_prepared_pdb()
    output:
        lowest_energy = f"{output_dir}/ligand_sampling/lowest_energy_conformer.pdb",
        conformer_list = f"{output_dir}/ligand_sampling/conformer_list.json"
    params:
        output_dir = f"{output_dir}/ligand_sampling",
        ligand_smiles = ligand_smiles,
        num_conformers = num_conformers,
        num_positions = num_positions,
        conformer_rmsd_cutoff = conformer_rmsd_cutoff,
        rotation = rotation,
        translation = translation,
        ligand_name = ligand_name
    run:
        from ligand_sampling import run_ligand_sampling_pipeline
        
        Path(params.output_dir).mkdir(parents=True, exist_ok=True)
        
        # Run ligand sampling
        result = run_ligand_sampling_pipeline(
            structure_path=input.pdb,
            ligand_smiles=params.ligand_smiles,
            output_path=params.output_dir,
            num_conformers=params.num_conformers,
            num_positions=params.num_positions,
            conformer_rmsd_cutoff=params.conformer_rmsd_cutoff,
            rotation=params.rotation,
            translation=params.translation,
            ligand_name=params.ligand_name
        )
        
        # Collect all generated conformer PDB files (num_conformers × num_positions)
        import glob
        import os
        import shutil
        
        # Check multiple possible locations for the output PDB files
        final_pdbs_dir = os.path.join(params.output_dir, 'final_pdbs')
        
        # Try to find PDB files in final_pdbs subdirectory first
        if os.path.isdir(final_pdbs_dir):
            conformer_pdbs = sorted(glob.glob(f"{final_pdbs_dir}/*.pdb"))
            print(f"  Found {len(conformer_pdbs)} PDB files in final_pdbs/")
        else:
            # Fallback: look in main directory
            conformer_pdbs = sorted(glob.glob(f"{params.output_dir}/*.pdb"))
            print(f"  Found {len(conformer_pdbs)} PDB files in main directory")
        
        # Remove the lowest energy file from the list if it exists (it will be separate)
        conformer_pdbs = [pdb for pdb in conformer_pdbs if 'lowest_energy' not in pdb]
        
        # Find the lowest energy conformer
        lowest_energy_pdb = None
        if os.path.isdir(final_pdbs_dir):
            lowest_energy_candidates = glob.glob(f"{final_pdbs_dir}/*lowest_energy*.pdb")
            if lowest_energy_candidates:
                lowest_energy_pdb = lowest_energy_candidates[0]
        
        # If no lowest energy found, use the first one if available
        if not lowest_energy_pdb and conformer_pdbs:
            lowest_energy_pdb = conformer_pdbs[0]
        
        # Copy or create lowest_energy_conformer.pdb in the main output directory
        if lowest_energy_pdb and os.path.exists(lowest_energy_pdb):
            shutil.copy(lowest_energy_pdb, output.lowest_energy)
            print(f"✓ Copied lowest energy conformer to {output.lowest_energy}")
        else:
            # Create a dummy file if no lowest energy found
            print(f"⚠ Warning: No lowest energy conformer found, creating placeholder")
            with open(output.lowest_energy, 'w') as f:
                f.write("# Placeholder - no lowest energy conformer found\n")
        
        # Save list of all conformers to JSON
        conformer_data = {
            'total_conformers': len(conformer_pdbs),
            'num_conformers': params.num_conformers,
            'num_positions': params.num_positions,
            'expected_total': params.num_conformers * params.num_positions,
            'conformers': conformer_pdbs,
            'lowest_energy': output.lowest_energy
        }
        
        with open(output.conformer_list, 'w') as f:
            json.dump(conformer_data, f, indent=2)
        
        print(f"✓ Ligand sampling completed. Results in {params.output_dir}")
        print(f"  Generated {len(conformer_pdbs)} PDB files ({params.num_conformers} conformers × {params.num_positions} positions)")
        if conformer_pdbs:
            print(f"  PDB files found: {[os.path.basename(p) for p in conformer_pdbs[:3]]}{'...' if len(conformer_pdbs) > 3 else ''}")
        print(f"  Conformer list saved to {output.conformer_list}")

# ============================================================================
# 4. GENERATE CONTIG MAP
# ============================================================================

rule generate_contig_map:
    """
    Generate the contig map for RF Diffusion.
    Defines which regions of the protein will be redesigned.
    """
    input:
        pdb_info = get_pdb_info_file()
    output:
        contig_map_file = f"{output_dir}/PDB_info/contig_map.json"
    params:
        output_dir = f"{output_dir}/PDB_info",
        segment_extension = segment_extension,
        n_termini_extension = n_termini_extension,
        c_termini_extension = c_termini_extension,
        conservative_redesign = conservative_redesign,
        user_defined_contig_map = user_defined_contig_map
    run:
        # Load PDB info
        with open(input.pdb_info, 'r') as f:
            pdb_info = json.load(f)
        
        # Generate contig map
        contig_map, segments = func.list_to_contig_map(
            chain_id=pdb_info['chain_id'],
            seq_length=int(pdb_info['sequence_length']),
            active_site=eval(pdb_info['active_site']),
            missing_residues=eval(pdb_info['missing_positions']),
            start=int(pdb_info['start_residue_number']),
            segment_extension=params.segment_extension,
            n_termini_extension=params.n_termini_extension,
            c_termini_extension=params.c_termini_extension,
            conservative_RF=params.conservative_redesign,
            DSSP_string=pdb_info['dssp_string'],
            user_defined_contig_map=params.user_defined_contig_map
        )
        
        # Save contig map and segments
        output_data = {
            'contig_map': contig_map,
            'segments': segments,
            'chain_id': pdb_info['chain_id'],
            'sequence_length': pdb_info['sequence_length']
        }
        
        with open(output.contig_map_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        print(f"✓ Contig map generated: {contig_map}")

# ============================================================================
# 5. RUN RF DIFFUSION
# ============================================================================

rule run_rf_diffusion:
    """
    Run RF Diffusion to generate binding site redesigns for each conformer.
    Supports three models: all_atom, RF1, or RF3.
    Processes all ligand conformers from the sampling step.
    """
    input:
        pdb = get_prepared_pdb(),
        contig_map_info = f"{output_dir}/PDB_info/contig_map.json",
        conformer_list = f"{output_dir}/ligand_sampling/conformer_list.json",
        pdb_info = get_pdb_info_file()
    output:
        designs_done = f"{output_dir}/RF_designs/designs_completed.txt"
    params:
        output_dir = f"{output_dir}/RF_designs",
        rf_model = RF_model,
        ligand_resname = ligand_resname,
        num_designs = num_designs,
        T = T,
        deterministic = deterministic,
        design_startnum = design_startnum,
        # All atom specific
        path_to_RFAA_apptainer = path_to_RFAA_apptainer,
        path_to_RFAA_script = path_to_RFAA_script,
        path_to_RFAA_weights = path_to_RFAA_weights,
        RFAA_ligand_name = RFAA_ligand_name,
        # RF1 specific
        path_to_RF1_script = path_to_RF1_script,
        path_to_RF1_env = path_to_RF1_env,
        # RF3 specific
        path_to_RF3_env = path_to_RF3_env,
        RF3_checkpoint_path = RF3_checkpoint_path,
        RF3_batch_size = RF3_batch_size,
        RF3_use_classifier_free_guidance = RF3_use_classifier_free_guidance,
        RF3_cfg_scale = RF3_cfg_scale,
        RF3_num_timesteps = RF3_num_timesteps,
        RF3_step_scale = RF3_step_scale,
        RF3_noise_scale = RF3_noise_scale,
        RF3_gamma_0 = RF3_gamma_0,
        RF3_gamma_min = RF3_gamma_min,
        RF3_dump_trajectories = RF3_dump_trajectories,
        RF3_prevalidate_inputs = RF3_prevalidate_inputs,
        RF3_low_memory_mode = RF3_low_memory_mode
    resources:
        gpus=1
    run:
        Path(params.output_dir).mkdir(parents=True, exist_ok=True)
        
        # Load contig map info and PDB info
        with open(input.contig_map_info, 'r') as f:
            contig_info = json.load(f)
        contig_map = contig_info['contig_map']
        chain_id = contig_info['chain_id']
        
        with open(input.pdb_info, 'r') as f:
            pdb_info = json.load(f)
        
        # Load list of all conformers to process
        with open(input.conformer_list, 'r') as f:
            conformer_data = json.load(f)
        
        # Extract conformer PDB paths - handle different possible JSON structures
        if isinstance(conformer_data, dict) and 'conformers' in conformer_data:
            conformer_pdbs = conformer_data['conformers']
        elif isinstance(conformer_data, list):
            conformer_pdbs = conformer_data
        else:
            # Fallback: try to get PDB paths from dict values
            conformer_pdbs = [v for v in conformer_data.values() if isinstance(v, str)]
        
        if not conformer_pdbs:
            raise ValueError(f"No conformer PDBs found in {input.conformer_list}")
        
        # Process each conformer
        completed_conformers = []
        for idx, conformer_pdb in enumerate(conformer_pdbs):
            # Create subdirectory for this conformer's designs
            conformer_output_dir = f"{params.output_dir}/conformer_{idx}"
            Path(conformer_output_dir).mkdir(parents=True, exist_ok=True)
            
            print(f"\n{'='*60}")
            print(f"Running RF Diffusion on conformer {idx} / {len(conformer_pdbs)}")
            print(f"Input PDB: {conformer_pdb}")
            print(f"Output: {conformer_output_dir}")
            print(f"{'='*60}\n")
            
            try:
                # Run appropriate RF diffusion model
                if params.rf_model == 'all_atom':
                    func.run_rfAA(
                        output_path=conformer_output_dir,
                        input_pdb=conformer_pdb,
                        contig_map=contig_map,
                        num_designs=params.num_designs,
                        T=params.T,
                        path_to_RFAA_apptainer=params.path_to_RFAA_apptainer,
                        path_to_RFAA_script=params.path_to_RFAA_script,
                        path_to_RFAA_weights=params.path_to_RFAA_weights,
                        inference_ligand=params.RFAA_ligand_name,
                        design_startnum=params.design_startnum,
                        deterministic=params.deterministic
                    )
                
                elif params.rf_model == 'RF1':
                    func.run_rf1(
                        output_path=conformer_output_dir,
                        input_pdb=conformer_pdb,
                        contig_map=contig_map,
                        num_designs=params.num_designs,
                        T=params.T,
                        path_to_RF1_script=params.path_to_RF1_script,
                        path_to_RF1_env=params.path_to_RF1_env
                    )
                
                elif params.rf_model == 'RF3':
                    func.run_rfd3(
                        output_path=conformer_output_dir,
                        input_pdb=conformer_pdb,
                        contig_map=contig_map,
                        pdb_info=pdb_info,
                        num_designs=params.num_designs,
                        chain_id=chain_id,
                        path_to_RF3_env=params.path_to_RF3_env,
                        ligand_resname=params.ligand_resname,
                        checkpoint_path=params.RF3_checkpoint_path,
                        batch_size=params.RF3_batch_size,
                        use_classifier_free_guidance=params.RF3_use_classifier_free_guidance,
                        cfg_scale=params.RF3_cfg_scale,
                        num_timesteps=params.RF3_num_timesteps,
                        step_scale=params.RF3_step_scale,
                        noise_scale=params.RF3_noise_scale,
                        gamma_0=params.RF3_gamma_0,
                        gamma_min=params.RF3_gamma_min,
                        dump_trajectories=params.RF3_dump_trajectories,
                        prevalidate_inputs=params.RF3_prevalidate_inputs,
                        low_memory_mode=params.RF3_low_memory_mode
                    )
                
                completed_conformers.append({
                    'conformer_idx': idx,
                    'conformer_pdb': conformer_pdb,
                    'output_dir': conformer_output_dir,
                    'status': 'completed'
                })
                print(f"✓ Conformer {idx} completed successfully")
            
            except Exception as e:
                print(f"✗ Conformer {idx} failed with error: {str(e)}")
                completed_conformers.append({
                    'conformer_idx': idx,
                    'conformer_pdb': conformer_pdb,
                    'output_dir': conformer_output_dir,
                    'status': 'failed',
                    'error': str(e)
                })
        
        # Create completion summary
        with open(output.designs_done, 'w') as f:
            f.write(f"RF Diffusion designs completed using {params.rf_model} model\n")
            f.write(f"Total conformers processed: {len(completed_conformers)}\n")
            f.write(f"Successfully completed: {sum(1 for c in completed_conformers if c['status'] == 'completed')}\n")
            f.write(f"Failed: {sum(1 for c in completed_conformers if c['status'] == 'failed')}\n")
            f.write(f"Number of designs per conformer: {params.num_designs}\n")
            if params.rf_model != 'RF3':
                f.write(f"Diffusion steps (T): {params.T}\n")
            f.write(f"\nDetailed results:\n")
            for result in completed_conformers:
                f.write(f"  - Conformer {result['conformer_idx']}: {result['status']} - {result['output_dir']}\n")
        
        print(f"\nRF Diffusion completed with {params.rf_model} model")
        print(f"Total: {len(completed_conformers)} conformers, "
              f"Completed: {sum(1 for c in completed_conformers if c['status'] == 'completed')}, "
              f"Failed: {sum(1 for c in completed_conformers if c['status'] == 'failed')}")


# ============================================================================
# 6. UNZIP RF DESIGN ARCHIVES INTO FINAL CIF FOLDER
# ============================================================================

rule unzip_rf_design_archives:
    """
    Find every .zip and .cif.gz file recursively under RF_designs,
    extract/decompress each archive into RF_final_cif, and record completion.
    """
    input:
        designs_done = f"{output_dir}/RF_designs/designs_completed.txt"
    output:
        unzip_done = f"{output_dir}/RF_final_cif/unzip_completed.txt"
    params:
        rf_designs_dir = f"{output_dir}/RF_designs",
        rf_final_cif_dir = f"{output_dir}/RF_final_cif"
    run:
        Path(params.rf_final_cif_dir).mkdir(parents=True, exist_ok=True)

        zip_files = []
        cif_gz_files = []
        for root, _, files in os.walk(params.rf_designs_dir):
            for file_name in files:
                if file_name.lower().endswith('.zip'):
                    zip_files.append(os.path.join(root, file_name))
                elif file_name.lower().endswith('.cif.gz'):
                    cif_gz_files.append(os.path.join(root, file_name))

        zip_files = sorted(zip_files)
        cif_gz_files = sorted(cif_gz_files)

        extracted_zip_count = 0
        decompressed_cif_gz_count = 0

        # Extract .zip archives
        for zip_path in zip_files:
            rel_parent = os.path.relpath(os.path.dirname(zip_path), params.rf_designs_dir)
            zip_stem = Path(zip_path).stem

            # Keep extracted outputs organized and avoid filename collisions.
            if rel_parent == '.':
                extract_dir = os.path.join(params.rf_final_cif_dir, zip_stem)
            else:
                extract_dir = os.path.join(params.rf_final_cif_dir, rel_parent, zip_stem)

            Path(extract_dir).mkdir(parents=True, exist_ok=True)
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                zip_ref.extractall(extract_dir)

            extracted_zip_count += 1
            print(f"✓ Extracted: {zip_path} -> {extract_dir}")

        # Decompress .cif.gz files into .cif files
        for cif_gz_path in cif_gz_files:
            rel_parent = os.path.relpath(os.path.dirname(cif_gz_path), params.rf_designs_dir)
            cif_name = Path(cif_gz_path).name[:-3]  # remove trailing '.gz'

            if rel_parent == '.':
                out_dir = params.rf_final_cif_dir
            else:
                out_dir = os.path.join(params.rf_final_cif_dir, rel_parent)

            Path(out_dir).mkdir(parents=True, exist_ok=True)
            out_cif_path = os.path.join(out_dir, cif_name)

            with gzip.open(cif_gz_path, 'rb') as src, open(out_cif_path, 'wb') as dst:
                shutil.copyfileobj(src, dst)

            decompressed_cif_gz_count += 1
            print(f"✓ Decompressed: {cif_gz_path} -> {out_cif_path}")

        with open(output.unzip_done, 'w') as f:
            f.write("RF archive extraction completed\n")
            f.write(f"Source directory: {params.rf_designs_dir}\n")
            f.write(f"Destination directory: {params.rf_final_cif_dir}\n")
            f.write(f"Zip files found: {len(zip_files)}\n")
            f.write(f"Zip files extracted: {extracted_zip_count}\n")
            f.write(f"CIF.GZ files found: {len(cif_gz_files)}\n")
            f.write(f"CIF.GZ files decompressed: {decompressed_cif_gz_count}\n")

        print(
            f"\nRF archive extraction complete. "
            f"Found {len(zip_files)} zip file(s) and {len(cif_gz_files)} cif.gz file(s)."
        )




