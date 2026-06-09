# Utility Scripts

This folder contains various utility scripts used throughout the project for tasks like data processing, PDB manipulation, structure analysis, and sequence optimization.

---

## Scripts overview

| Script | Purpose |
|--------|---------|
| [`calculate_sap_score.py`](#calculate_sap_scorepy) | Calculate SAP aggregation scores for PDB structures |
| [`conda_exporter.sh`](#conda_exportersh) | Export all local Conda environments to `.yml` files |
| [`extract_monomer.py`](#extract_monomerpy) | Extract Chain A monomer from multimeric PDB/CIF structures |
| [`extract_seq_pdb.py`](#extract_seq_pdbpy) | Extract amino acid sequences from PDB files into a CSV |
| [`fasta2xlsx.py`](#fasta2xlsxpy) | Convert a FASTA file into a 96-well plate Excel layout |
| [`hth_family_screening.py`](#hth_family_screeningpy) | Screen RCSB PDB for LacI-family HTH candidates |
| [`ligand_sampling.py`](#ligand_samplingpy) | Sample novel ligand conformations into a reference binding site |
| [`protein_to_codon_opt_DNA.py`](#protein_to_codon_opt_dnapy) | Generate codon-optimized DNA sequences via IDT's API |

---

## Script details

### `calculate_sap_score.py`
Calculates the Spatial Aggregation Propensity (SAP) score for a directory of PDB structures. It evaluates the hydrophobicity and solvent accessibility of residues within a default radius of 5.0 Å to score potential aggregation-prone regions, and outputs results to a CSV file.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input_folder` | Path to the folder containing PDB files |
| `--output_folder` | Path to the folder where `sap_results.csv` will be saved |

```bash
python calculate_sap_score.py --input_folder path/to/pdbs --output_folder path/to/output
```

---

### `conda_exporter.sh`
A bash utility that iterates through all local Conda environments and exports them to `.yml` files using `--no-builds` for cross-platform compatibility.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `TARGET_DIR` (positional) | Destination directory for exported environments (defaults to `./conda_exports`) |

```bash
./conda_exporter.sh /path/to/destination
```

---

### `extract_monomer.py`
Uses PyMOL to parse multimeric protein structures (`.pdb` or `.cif`) and extracts only Chain A, saving it as a new PDB file.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input_folder` | Path to the folder containing multimeric PDB/CIF files |
| `--output_folder` | Path to the folder where monomer PDB files will be saved |

```bash
python extract_monomer.py --input_folder path/to/multimers --output_folder path/to/monomers
```

---

### `extract_seq_pdb.py`
Reads all PDB files in a directory, extracts their amino acid sequences using BioPython's Polypeptide builder, and compiles them into a CSV mapping protein name to sequence.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input_folder` | Path to the folder containing PDB files |
| `--output_file` | Path and filename for the output CSV |

```bash
python extract_seq_pdb.py --input_folder path/to/pdbs --output_file results.csv
```

---

### `fasta2xlsx.py`
Parses a FASTA file, performs padding and standard N/C termini additions (e.g. to reach a minimum 300 nt length), and maps sequences into a 96-well plate layout exported as an Excel `.xlsx` file.

> **Note:** Input/output paths are currently hardcoded in the `__main__` block. Edit the script to point to your specific files before running.

```bash
python fasta2xlsx.py
```

---

### `hth_family_screening.py`
Mines the RCSB PDB via its GraphQL API to identify all structures belonging to the LacI family (PF00356). Groups them by UniProt accession and analyzes pairs of structures (e.g. apo vs. holo states) to measure conformational hinge shifts, helping identify suitable chassis candidates for de novo allostery design.

> No arguments required. Parameters and target families are configured internally.

```bash
python hth_family_screening.py
```

---

### `ligand_sampling.py`
A pre-processing step for binding site redesign. Given a novel ligand (as a SMILES string) and a reference protein structure with an existing ligand bound, it uses RDKit to sample new conformers and PyMOL to sample local positions and rotations. Outputs PDB files of the substituted complexes.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--ligand` | SMILES string of the novel ligand |
| `--structure` | Path to the reference PDB/CIF structure containing the original ligand |
| `--output` | Directory where output conformations and complexes will be saved |
| `--n_confs` *(optional)* | Number of RDKit conformers to generate (default: `5`) |
| `--n_pos` *(optional)* | Number of positional variations to sample (default: `5`) |

```bash
python ligand_sampling.py --ligand "CC(=O)OC1=CC=CC=C1C(=O)O" --structure complex.pdb --output ./sampled_complexes
```

---

### `protein_to_codon_opt_DNA.py`
Takes a CSV of protein sequences and uses IDT's API to generate codon-optimized DNA sequences for *Escherichia coli* K12. Handles restriction site avoidance, limits cysteine count, adjusts net charge, and adds required flanking sequences. Outputs a multi-FASTA file for both processed proteins and optimized DNA.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input_csv` | Path to the input CSV (must contain `Protein Name` and `Sequence` columns) |
| `--output_dir` | Directory to save the resulting CSV and `.fasta` files |

```bash
python protein_to_codon_opt_DNA.py --input_csv sequences.csv --output_dir ./optimized_dna
```
