# RF Diffusion Binding Site Designer - Snakemake Pipeline

## Overview

This Snakemake pipeline implements a binding site redesign workflow using RF Diffusion and LigandMPNN. The pipeline automates the process of:

1. **PDB Preparation**: Cleaning and preparing PDB structures
2. **PDB Analysis**: Extracting sequence, secondary structure, and defining active sites
3. **Ligand Sampling**: Generating ligand conformers and exploring binding pocket positions
4. **Contig Map Generation**: Defining protein regions for redesign
5. **RF Diffusion**: Generating binding site redesigns (supports all_atom, RF1, or RF3 models)
6. **LigandMPNN Integration**: (Optional) Sequence design for redesigned structures

## Pipeline Architecture

```
structure_path
    ↓
[clean_pdb]  ← Remove water/non-protein atoms
    ↓
    ├─→ [extract_pdb_info] ← Extract sequence, DSSP, active site
    │       ↓
    │   pdb_info.json
    │       ↓
    │   [generate_contig_map] ← Define RF diffusion regions
    │       ↓
    │   contig_map.json
    │
    └─→ [ligand_sampling] ← Generate conformers & positions
            ↓
        ligand_conformers.pdb
        ↓
        [run_rf_diffusion] ← Run RF diffusion, choose model
        ↓
    RF_designs/ (multiple PDB designs)
```

## Installation & Setup

### Requirements

- **Python 3.8+**
- **Snakemake 7.0+**
- **Conda/Mamba** (for environment management)
- **Required packages**:
  - Biopython
  - NumPy, Pandas
  - PyMOL (for ligand detection)
  - RDKit (for ligand handling)
  - tmtools (for TM-align)

### External Tools

- **RF Diffusion All Atom** (optional - for all_atom model)
- **RF1 or RF3** (optional - if using those models)
- **LigandMPNN** (optional - for sequence design)
- **ESMFold** (optional - for structure validation)

### Installation Steps

```bash
# 1. Clone or navigate to the pipeline directory
cd allosteric_biosensor_design_prediction/strategy_1_snakemake

# 2. Create conda environment (optional)
conda create -n rf_diffusion_bsd -c conda-forge snakemake python=3.10
conda activate rf_diffusion_bsd

# 3. Install dependencies
pip install -r requirements.txt  # if available
# OR manually:
pip install biopython numpy pandas rdkit pymol tmtools

# 4. Update configuration file (see below)
nano config_rfdiffusion_bsd.yml
```

## Configuration

### Main Configuration File: `config_rfdiffusion_bsd.yml`

The pipeline is controlled by a YAML configuration file with the following sections:

#### 1. **Basic Parameters**
```yaml
structure_path: "path/to/input.pdb"          # Input structure
ligand_smiles: "CCO"                          # Target ligand SMILES
ligand_name: "Tc"                             # Ligand name in PDB
```

#### 2. **Ligand Sampling Parameters**
```yaml
num_conformers: 2                             # Number of conformers to generate
num_positions: 1                              # Number of sampling positions
conformer_rmsd_cutoff: 0.75                   # Conformer diversity threshold (Å)
rotation: 15                                  # Rotation angle (degrees)
translation: 0.3                              # Translation distance (Å)
```

#### 3. **Active Site Definition**
```yaml
first_shell_distance: 5.0                     # Distance to define first shell (Å)
second_shell_distance: 5.0                    # Distance to define second shell (Å)
include_second_shell: false                   # Include second shell in redesign
user_defined_active_site: null                # Override with [10, 15, 20]
user_defined_residues: null                   # Residues to force include [5, 8]
```

#### 4. **Redesign Conditions**
```yaml
RF_model: "all_atom"                          # Model: "all_atom", "RF1", or "RF3"
num_designs: 10                               # Number of designs to generate
T: 50                                         # Diffusion steps
conservative_redesign: false                  # Conserve secondary structure
segment_extension: 1                          # N/C extension from active site (residues)
n_termini_extension: 1                        # N-terminus extension
c_termini_extension: null                     # C-terminus extension
user_defined_contig_map: null                 # Override contig map
```

#### 5. **RF Diffusion Parameters** (if using all_atom model)
```yaml
path_to_RFAA_apptainer: "/path/to/rf_se3_diffusion.sif"
path_to_RFAA_script: "/path/to/run_inference.py"
path_to_RFAA_weights: "/path/to/RFDiffusionAA_paper_weights.pt"
RFAA_ligand_name: "UNL"                       # Ligand name in RF model
deterministic: true                           # Reproducible results
design_startnum: 1                            # Starting design number
```

#### 6. **LigandMPNN Parameters** (for sequence design)
```yaml
path_to_ligand_MPNN_script: "/path/to/LigandMPNN/run.py"
path_to_ligand_MPNN_env: "/path/to/ligandmpnn_env"
MPNN_num_designs: 10                          # Sequences per structure
mpnn_temperature: 0.05                        # Sampling temperature
first_shell_only: false                       # Only redesign first shell
```

## Usage

### Basic Execution

```bash
# Dry run (check DAG without executing)
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -n

# Run pipeline with 1 job
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1

# Run with 4 parallel jobs
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 4

# Run with verbose output
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 -p

# Generate DAG visualization
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  --dag | dot -Tpng > dag.png
```

### Advanced Options

```bash
# Run specific rule
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 run_rf_diffusion

# Print shell commands for debugging
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 --printshellcmds

# Remove specific outputs to re-run
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 --forcerun run_rf_diffusion

# Use conda environments (if defined)
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 --use-conda
```

## Output Structure

```
output_dir/
├── PDB_prep/
│   ├── cleaned_structure.pdb          # Cleaned input structure
│   └── clean_status.txt               # Status file
├── PDB_info/
│   ├── pdb_info.json                  # Extracted PDB information
│   └── contig_map.json                # RF diffusion contig map
├── ligand_sampling/
│   ├── lowest_energy_conformer.pdb    # Best ligand conformer
│   └── conformer_list.json            # All generated conformers
└── RF_designs/
    ├── design_*.pdb                   # Individual RF diffusion designs
    ├── scores.json                    # Design evaluation metrics
    └── designs_completed.txt           # Completion marker
```

## Pipeline Rules

### Rule: `clean_pdb`
- **Input**: Original PDB file
- **Output**: Cleaned PDB (protein + ligand only)
- **Function**: Removes water, ions, and other non-relevant atoms

### Rule: `extract_pdb_info`
- **Input**: Cleaned PDB
- **Output**: JSON with sequence, DSSP, active site information
- **Function**: Uses `func.extract_pdb_info()` from functions_bsd.py

### Rule: `ligand_sampling`
- **Input**: Cleaned PDB, ligand SMILES
- **Output**: Sampled ligand conformers and positions
- **Function**: Runs `run_ligand_sampling_pipeline()` from ligand_sampling.py

### Rule: `generate_contig_map`
- **Input**: PDB info JSON
- **Output**: Contig map JSON for RF diffusion
- **Function**: Uses `func.list_to_contig_map()` to define redesign regions

### Rule: `run_rf_diffusion`
- **Input**: Prepared PDB, contig map, ligand conformer
- **Output**: RF diffusion designs
- **Function**: Calls appropriate RF model (`run_rfAA()`, `run_rf1()`, `run_rf3()`)
- **Supported Models**:
  - `all_atom`: RF Diffusion All Atom (recommended for small molecules)
  - `RF1`: RF Diffusion 1 (protein-only)
  - `RF3`: RF Diffusion 3 (next-gen model)

## Workflow Notes

### Important Considerations

1. **Ligand Position**: Ensure the input PDB has the original ligand in the binding pocket for active site detection
2. **GPU Requirements**: RF Diffusion typically requires GPU acceleration (NVIDIA with CUDA)
3. **Path Configuration**: Update all paths in config to match your system setup
4. **Deterministic Mode**: Set `deterministic: true` for reproducible results
5. **Design Number**: Each design is independently generated; increase `num_designs` for more options

### Extending the Pipeline

To add LigandMPNN integration after RF Diffusion:

1. Uncomment the `generate_mpnn_jsons` and `run_ligand_mpnn` rules (at end of snakemake file)
2. Implement the function calls using `func.generate_MPNN_jsons()` and `func.run_ligand_MPNN()`
3. Update `rule all` to include MPNN outputs
4. Add required MPNN parameters to config file

## Troubleshooting

### Common Issues

| Issue | Solution |
|-------|----------|
| `FileNotFoundError: structure_path` | Check path in config file exists |
| `ModuleNotFoundError: functions_bsd` | Ensure functions_bsd.py is in same directory |
| `DSSP not found` | Install DSSP: `sudo apt install dssp` |
| `PyMOL command failed` | Install PyMOL with conda: `conda install pymol` |
| `RF Diffusion not found` | Update paths in config to correct locations |
| `CUDA out of memory` | Reduce `num_designs` or run on GPU with more memory |

### Debugging

```bash
# Run with debug output
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 --debug

# Check log files in output directory
tail -f outputs/rf_diffusion_bsd/RF_designs/*.log

# Validate config syntax
python -c "import yaml; print(yaml.safe_load(open('config_rfdiffusion_bsd.yml')))"
```

## References

- **RF Diffusion**: https://github.com/RosettaCommons/RFdiffusion
- **LigandMPNN**: https://github.com/ProteinMPNN/LigandMPNN
- **Snakemake Documentation**: https://snakemake.readthedocs.io/
- **Biopython**: https://biopython.org/

## Citation

If you use this pipeline, please cite:

```
RF Diffusion paper (insert citation)
LigandMPNN paper (insert citation)
Snakemake paper (if applicable)
```

## Contact & Support

For issues, questions, or contributions, please contact:
- Eduardo Amo González (original author of binding_site_designer.py)
- [Your contact information]

---

**Last Updated**: March 2026
