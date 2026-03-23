# RF Diffusion Binding Site Designer - Quick Start Guide

## 5-Minute Setup

### Step 1: Configure Your Input
Edit `config_rfdiffusion_bsd.yml`:

```yaml
# REQUIRED - Update these three lines
structure_path: "/path/to/your/protein_with_ligand.pdb"
ligand_smiles: "CC(=O)O"  # Your target ligand SMILES
ligand_name: "ACE"         # 3-letter code for your ligand

# RECOMMENDED - Update these paths
path_to_RFAA_apptainer: "/path/to/rf_se3_diffusion.sif"
path_to_RFAA_script: "/path/to/run_inference.py"
path_to_RFAA_weights: "/path/to/RFDiffusionAA_paper_weights.pt"

# Optional - Adjust defaults based on your protein/ligand
num_designs: 10            # Number of designs (10-100)
T: 50                      # Diffusion steps (50-200)
first_shell_distance: 5.0  # Define active site (Å)
```

### Step 2: Verify Installation

```bash
# Check snakemake
python -c "import snakemake; print('✓ Snakemake ready')"

# Check dependencies
python -c "import Bio, numpy, pandas, rdkit; print('✓ Dependencies ready')"

# Verify PyMOL
python -c "from pymol import cmd; print('✓ PyMOL ready')"
```

### Step 3: Run Dry Run

```bash
cd strategy_1_snakemake

snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -n --quiet
```

**Expected Output:**
```
Job counts:
  clean_pdb: 1
  extract_pdb_info: 1
  ligand_sampling: 1
  generate_contig_map: 1
  run_rf_diffusion: 1
  all: 1
```

If you see this, you're ready to run!

### Step 4: Run Full Pipeline

```bash
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  -c 1 -p
```

Optional: Use `-c 4` for 4 parallel jobs (if your system supports it)

## What Each Rule Does

| Rule | Input | Output | Time |
|------|-------|--------|------|
| `clean_pdb` | PDB file | Cleaned PDB | <1 min |
| `extract_pdb_info` | Cleaned PDB | JSON with sequence info | <1 min |
| `ligand_sampling` | Cleaned PDB + SMILES | Ligand conformers | 1-5 min |
| `generate_contig_map` | PDB info | Redesign regions | <1 min |
| `run_rf_diffusion` | All above | Designed structures | 10-60 min* |

*Depends on num_designs and T values

## Output Files

After successful run, check:

```
outputs/rf_diffusion_bsd/
├── PDB_prep/cleaned_structure.pdb      # Your cleaned input
├── PDB_info/pdb_info.json              # Analysis results
├── ligand_sampling/                    # Ligand data
└── RF_designs/                         # Final designs!
    ├── design_1.pdb
    ├── design_2.pdb
    └── ...
```

## Common Configurations

### Quick Test (2 min)
```yaml
num_designs: 1
T: 5
num_conformers: 1
num_positions: 1
```

### Fast Run (10 min)
```yaml
num_designs: 5
T: 25
num_conformers: 2
num_positions: 1
```

### Production Run (1 hour)
```yaml
num_designs: 50
T: 100
num_conformers: 5
num_positions: 3
conservative_redesign: true
```

## Troubleshooting

### Error: "structure_path not found"
```bash
# Check your path is correct
ls /path/to/your/protein.pdb
# Update config_rfdiffusion_bsd.yml with correct path
```

### Error: "ModuleNotFoundError: functions_bsd"
```bash
# Ensure you're running from the correct directory
cd strategy_1_snakemake
pwd  # Should end with strategy_1_snakemake
```

### Error: "RF Diffusion not found"
```bash
# Check paths in config
ls /path/to/rf_se3_diffusion.sif
ls /path/to/run_inference.py
# Update if paths are incorrect
```

### Long Waits (typical)
- First run of `clean_pdb`: Normal, reading large structure
- `ligand_sampling`: 1-5 min depending on conformers
- `run_rf_diffusion`: Expected! This is the main computation (10-60 min)

## Next Steps

### View Pipeline DAG
```bash
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  --dag | dot -Tpng > pipeline.png
```

### Run Specific Rule
```bash
# Rerun just the PDB info extraction
snakemake --snakefile rf_diffusion_bsd_snakemake.smk \
  --configfile config_rfdiffusion_bsd.yml \
  extract_pdb_info -c 1 -p
```

### Add LigandMPNN Integration
After RF Diffusion completes:

1. Uncomment MPNN rules at end of `rf_diffusion_bsd_snakemake.smk`
2. Update MPNN paths in config
3. Rerun pipeline: Snakemake will only run new rules

### Analyze Results
```python
import json
import pandas as pd

# Load design info
with open('outputs/rf_diffusion_bsd/RF_designs/scores.json') as f:
    scores = json.load(f)

df = pd.DataFrame(scores)
print(df.sort_values('score'))
```

## Performance Tips

### Speed Up
- Reduce `num_designs` (fewer designs = faster)
- Reduce `T` (fewer diffusion steps = faster, lower quality)
- Use `first_shell_only: true` to redesign smaller region

### Quality Up
- Increase `num_designs` (more options to choose from)
- Increase `T` (more diffusion steps = higher quality)
- Set `conservative_redesign: true`
- Increase `num_conformers` for better ligand sampling

### GPU Optimization
```yaml
# For A100 or H100 GPUs
num_designs: 100
T: 200
num_conformers: 10

# For smaller GPUs (RTX 3090)
num_designs: 20
T: 50
num_conformers: 3
```

## Support

For detailed documentation, see [README.md](README.md)

For issues with specific functions, check [functions_bsd.py](functions_bsd.py)

---

**Ready?** Run: `snakemake --snakefile rf_diffusion_bsd_snakemake.smk --configfile config_rfdiffusion_bsd.yml -c 1`

