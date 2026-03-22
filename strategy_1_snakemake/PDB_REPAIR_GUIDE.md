# PDB Repair Integration Guide

## Overview

The `repair_pdb.py` module validates and repairs PDB structures before passing them to LigandMPNN. It:

1. **Detects chain breaks** - Identifies gaps in residue numbering and CA-CA distances
2. **Fills missing residues** - Inserts Glycine residues to bridge gaps
3. **Updates redesign lists** - Includes filled residues in the LigandMPNN redesign set
4. **Generates reports** - Documents all changes for transparency

## Why This Is Needed

The error you encountered:
```
ValueError: shape mismatch: value array of shape (1,2328) could not be broadcast to indexing result of shape (2324,)
```

This happens when:
- RF Diffusion outputs have truncated chains or missing loops
- ProDy expects continuous residue numbering
- LigandMPNN generates sequences that don't match the PDB backbone length

## Integration into Snakemake

### Option 1: Add as a Pre-Processing Rule (Recommended)

Add this to `bsd_def.smk` **before** `rule generate_mpnn_jsons`:

```smakemake
# --- PRE-PROCESSING: Validate and repair PDBs ---
rule repair_pdb_inputs:
    input:
        pdb_paths_json = f"{OUTPUT_DIR}/MPNN_jsons/pdb_paths_multi.json",
        redesign_json = f"{OUTPUT_DIR}/MPNN_jsons/redesigned_residues_multi.json"
    output:
        repaired_pdb_json = f"{OUTPUT_DIR}/MPNN_jsons/pdb_paths_multi_repaired.json",
        repaired_redesign_json = f"{OUTPUT_DIR}/MPNN_jsons/redesigned_residues_multi_repaired.json",
        repair_report = f"{OUTPUT_DIR}/repair_summary.txt"
    params:
        repair_script = "repair_pdb.py"
    run:
        from repair_pdb import batch_repair_and_update_jsons
        batch_repair_and_update_jsons(
            input.pdb_paths_json,
            input.redesign_json,
            Path(output.repaired_pdb_json).parent
        )
```

Then update the `split_mpnn_jsons` rule to use the repaired JSONs:

```smakemake
rule split_mpnn_jsons:
    input:
        paths_json = f"{OUTPUT_DIR}/MPNN_jsons/pdb_paths_multi_repaired.json",  # CHANGED
        res_json = f"{OUTPUT_DIR}/MPNN_jsons/redesigned_residues_multi_repaired.json"  # CHANGED
    output:
        ...
```

### Option 2: Standalone Batch Processing

Run before the Snakemake pipeline:

```bash
python repair_pdb.py --batch \
  strategy_1_campaign/sbatch_trial/MPNN_jsons/pdb_paths_multi.json \
  strategy_1_campaign/sbatch_trial/MPNN_jsons/redesigned_residues_multi.json \
  strategy_1_campaign/sbatch_trial/pdb_repairs/
```

This creates:
- `pdb_paths_multi_repaired.json` - Updated PDB paths
- `redesigned_residues_multi_repaired.json` - Updated redesign lists with filled residues
- `repair_summary.txt` - Full repair report

Then update your pipeline to use the repaired JSONs.

## Usage Examples

### Single PDB Repair

```python
from repair_pdb import repair_pdb_and_update_redesign_list

result = repair_pdb_and_update_redesign_list(
    pdb_path='/path/to/structure.pdb',
    redesign_residues_dict={'A': [100, 105, 110]},
    output_dir='./repaired/'
)

print(result['report'])
# Outputs:
# - Original residues: 324
# - Filled residues: 4  
# - New total: 328
# - Filled positions: [('A', 156), ('A', 157), ('B', 45), ('B', 46)]
```

### Batch Repair (Command Line)

```bash
python repair_pdb.py --batch pdb_paths.json redesign_residues.json output_dir/
```

Output files:
- `pdb_paths_multi_repaired.json` - Updated file
- `redesigned_residues_multi_repaired.json` - Updated file with new residues
- `repair_summary.txt` - Detailed report

### Batch Repair (Python)

```python
from repair_pdb import batch_repair_and_update_jsons

results = batch_repair_and_update_jsons(
    pdb_paths_json='pdb_paths_multi.json',
    redesigned_residues_json='redesigned_residues_multi.json',
    output_dir='./repaired/'
)

for pdb_name, result in results.items():
    print(f"{pdb_name}: {result['filled_residues']} residues filled")
```

## How It Works

### 1. Chain Break Detection

Uses two criteria:

**Spatial Gap:** CA-CA distance > 4.0Å (default)
```python
distance = ||CA_n+1 - CA_n|| > 4.0 Å → gap
```

**Numbering Gap:** Residue number jump > 1
```python
resnum_n+1 - resnum_n > 1 → missing residues
```

### 2. Residue Filling

For each gap `[resnum_start, resnum_end]`:
- Creates Glycine residues (G = most flexible)
- Interpolates backbone coordinates between anchors
- High B-factors (100.0) mark filled residues

### 3. Redesign List Update

```python
# Original redesign list for structure A
{'A': [100, 110, 120]}

# After repair finds gaps at 115-116
# Updated list:
{'A': [100, 110, 115, 116, 120]}
```

LigandMPNN will now redesign these filled residues!

## Configuration Parameters

In `repair_pdb.py`, adjust detection thresholds:

```python
# Gap detection parameters
gap_threshold_distance = 4.0      # Angstroms
gap_threshold_resnum = 1          # missing residue jump
```

## Output Files

### pdb_paths_multi_repaired.json
```json
{
    "structure_1": "/path/to/structure_1_repaired.pdb",
    "structure_2": "/path/to/structure_2_repaired.pdb"
}
```

### redesigned_residues_multi_repaired.json
```json
{
    "structure_1": {
        "A": [100, 105, 110, 115, 116],
        "B": [50, 55]
    },
    "structure_2": {
        "A": [20, 25, 30]
    }
}
```

### repair_summary.txt
```
PDB REPAIR SUMMARY
================================================================================

Processed: 5 PDBs
Successful: 5
Failed: 0

structure_1:
  Success: True
  Original residues: 324
  Filled residues: 4
  Filled at: [('A', 156), ('A', 157), ('B', 45), ('B', 46)]

structure_2:
  Success: True
  Original residues: 250
  Filled residues: 0

...
```

## Troubleshooting

### Error: "Could not parse PDB"
- Ensure valid PDB format
- Check file exists and is readable

### Error: "No CA atoms found"
- PDB may have no backbone
- Verify structure integrity

### Filled residues not being redesigned
- Check that `updated_redesign_dict` is used in MPNN call
- Verify redesign list format matches expectations

## Next Steps

1. **Test on your data:**
   ```bash
   python repair_pdb.py --batch pdb_paths_multi.json redesigned_residues_multi.json test_output/
   ```

2. **Review repair_summary.txt** to see what was changed

3. **Integrate into workflow** using one of the options above

4. **Run MPNN** using updated JSONs

## Performance Notes

- Single PDB: ~0.1-0.5s
- 100 PDBs: ~10-50s (depends on structure size and gap count)
- Memory: Minimal (~100MB for typical structures)
