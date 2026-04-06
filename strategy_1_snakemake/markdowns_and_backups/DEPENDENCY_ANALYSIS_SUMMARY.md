# BSD Function Library - Comprehensive Recursive Dependency Analysis

**Analysis Date:** 2026-03-31  
**Source File:** `functions_snakemake_bsd.py`  
**Pipeline File:** `bsd_def.smk`  
**Total Functions:** 21  
**Total Used:** 21 (100%)  
**Total Unused:** 0

---

## Executive Summary

All 21 functions defined in the BSD function library are actively used by the Snakemake pipeline. The library consists of:
- **18 Root Functions** - directly called by `bsd_def.smk`
- **3 Helper Functions** - called transitively through root functions

### Root Functions (18)
These functions are directly invoked by the pipeline:

1. `boltz_yaml_generator` - Generates Boltz-compatible YAML configuration files
2. `check_cofold_validity` - Validates cofold predictions by checking ligand placement
3. `filter_dataframe_1D` - Removes sequences with repeating amino acids
4. `generate_redesign_string` - Identifies residues requiring redesign in RF diffusion structures
5. `plot_comparative_distributions` - Creates KDE distribution plots across multiple datasets
6. `plot_four_scatters_with_region` - Generates 2x2 grid of scatter plots with target regions
7. `process_AF3_folder` - Consolidates AF3 prediction outputs
8. `process_Boltz_folder` - Processes Boltz predictions and aggregates metrics
9. `process_MPNN_folder` - Processes LigandMPNN output folders
10. `process_chai_folder` - Consolidates Chai prediction outputs
11. `run_AlphaFold3` - Executes AlphaFold3 predictions
12. `run_chai` - Runs Chai protein folding (basic)
13. `run_chai_cofold` - Runs Chai cofold predictions (protein-ligand)
14. `run_chai_cofold_w_MSA` - Runs Chai cofold with MSA
15. `run_chai_w_MSA` - Runs Chai with MSA
16. `threed_filter_1_df` - First-pass 3D structure filtering
17. `threed_filter_2_df` - Second-pass multi-replicate 3D filtering
18. `threed_params_1_df` - Calculates 3D metrics (pLDDT, RMSD, TMscore, GNINA scores)

### Helper Functions (3)
These functions are called indirectly through root functions:

1. `process_MPNN_file` - Called by `process_MPNN_folder`
   - Processes individual MPNN output files
   - Parses FASTA headers to extract confidence scores
   - Filters top sequences

2. `get_centroid_distance` - Called by `check_cofold_validity`
   - Calculates centroid distance between ligand and active site
   - Returns validity flag based on distance cutoff

3. `global_score` - Called by `threed_filter_1_df` and `threed_filter_2_df`
   - Normalizes 3D metrics (pLDDT, TMscore, clashes)
   - Calculates weighted global score for structure selection

---

## Dependency Tree

```
boltz_yaml_generator
├─ (no internal dependencies)

check_cofold_validity
├─ get_centroid_distance ◄── HELPER FUNCTION

filter_dataframe_1D
├─ sliding_window_1D_filter ◄── EXTERNAL (from other module)

generate_redesign_string
├─ load_structure ◄── EXTERNAL (from biotite)

plot_comparative_distributions
├─ (no internal dependencies)

plot_four_scatters_with_region
├─ count_designs_in_region ◄── EXTERNAL

process_AF3_folder
├─ (no internal dependencies)

process_Boltz_folder
├─ (no internal dependencies)

process_MPNN_folder
├─ process_MPNN_file ◄── HELPER FUNCTION

process_chai_folder
├─ (no internal dependencies)

run_AlphaFold3
├─ AF_json_generator_w_MSA ◄── EXTERNAL
├─ AF_json_generator_cofold_w_MSA ◄── EXTERNAL
├─ AF_json_generator_wo_MSA ◄── EXTERNAL
├─ AF_json_generator_cofold ◄── EXTERNAL

run_chai
├─ (no internal dependencies)

run_chai_cofold
├─ (no internal dependencies)

run_chai_cofold_w_MSA
├─ (no internal dependencies)

run_chai_w_MSA
├─ (no internal dependencies)

threed_filter_1_df
├─ global_score ◄── HELPER FUNCTION

threed_filter_2_df
├─ global_score ◄── HELPER FUNCTION

threed_params_1_df
├─ extract_plddt ◄── EXTERNAL
├─ batch_tm_align ◄── EXTERNAL
├─ detect_clashes ◄── EXTERNAL
├─ extract_sequence ◄── EXTERNAL
├─ gnina_minimize_defined_box ◄── EXTERNAL
```

---

## Notes on External Dependencies

The following functions are referenced but defined in external modules:

1. **Structural Analysis:**
   - `extract_plddt` - Extracts pLDDT confidence scores
   - `batch_tm_align` - Calculates TMscore alignment metrics
   - `detect_clashes` - Detects steric clashes in structures
   - `extract_sequence` - Extracts protein sequences from structures

2. **Docking & Scoring:**
   - `gnina_minimize_defined_box` - Runs GNINA molecular docking

3. **Filtering & Metrics:**
   - `sliding_window_1D_filter` - Identifies repeating amino acid patterns
   - `count_designs_in_region` - Counts designs within target region

4. **AF3 Integration:**
   - `AF_json_generator_w_MSA`
   - `AF_json_generator_cofold_w_MSA`
   - `AF_json_generator_wo_MSA`
   - `AF_json_generator_cofold`

5. **Structure Loading:**
   - `load_structure` - From biotite library

---

## Function Statistics

| Metric | Value |
|--------|-------|
| Total Functions Defined | 21 |
| Root Functions (Direct Pipeline) | 18 |
| Helper Functions (Indirect) | 3 |
| Unused Functions | 0 |
| Code Coverage | 100% |
| Average Function Size | ~52 lines |

---

## Pipeline Integration Map

```
bsd_def.smk (Snakemake Workflow)
    │
    ├─→ Rule: generate_mpnn_jsons
    │   └─→ generate_redesign_string
    │
    ├─→ Rule: process_mpnn_output_and_filter
    │   ├─→ process_MPNN_folder → process_MPNN_file
    │   └─→ filter_dataframe_1D
    │
    ├─→ Rule: process_esm_batch
    │   └─→ threed_params_1_df
    │
    ├─→ Rule: first_3d_filter (checkpoint)
    │   └─→ threed_filter_1_df → global_score
    │
    ├─→ Rule: prep_boltz_yaml
    │   └─→ boltz_yaml_generator
    │
    ├─→ Rule: process_boltz_folder
    │   └─→ process_Boltz_folder
    │
    ├─→ Rule: run_second_prediction
    │   ├─→ run_chai
    │   ├─→ run_chai_w_MSA
    │   ├─→ run_chai_cofold
    │   ├─→ run_chai_cofold_w_MSA
    │   └─→ run_AlphaFold3
    │
    ├─→ Rule: process_second_prediction_folder
    │   ├─→ process_chai_folder
    │   └─→ process_AF3_folder
    │
    ├─→ Rule: gather_final_candidates
    │   ├─→ threed_params_1_df (3x)
    │   ├─→ check_cofold_validity → get_centroid_distance
    │   └─→ threed_filter_2_df → global_score
    │
    ├─→ Rule: plot_plddt_vs_tmscore
    │   └─→ plot_four_scatters_with_region
    │
    ├─→ Rule: plot_gnina
    │   └─→ plot_four_scatters_with_region
    │
    └─→ Rule: plot_metric_distributions
        └─→ plot_comparative_distributions
```

---

## Key Findings

1. **No Dead Code:** All 21 functions in the library are actively used
2. **Well-Structured:** Clear separation between:
   - Data processing functions (MPNN, 3D metrics)
   - Prediction execution functions (Chai, AF3, Boltz)
   - Filtering functions (1D, 3D single-pass, 3D multi-pass)
   - Visualization functions (scatter plots, distributions)
3. **Modular Design:** Helper functions (`global_score`, `get_centroid_distance`, `process_MPNN_file`) are reused across multiple pipeline stages
4. **External Dependencies:** Heavy reliance on external module functions for structural analysis (extract_plddt, detect_clashes, etc.)

---

## Maintenance Recommendations

1. **Refactoring Candidates:**
   - Consider extracting common 3D analysis patterns used in `threed_params_1_df` and `threed_filter_1_df`

2. **Testing Priorities:**
   - Core functions: `threed_filter_1_df`, `threed_filter_2_df` (multi-stage filtering logic)
   - Data processors: `process_MPNN_file`, `process_Boltz_folder` (data integrity)

3. **Documentation:**
   - Add docstrings to helper functions emphasizing their role in the pipeline
   - Document external dependencies and their expected behavior

