# Scripts

This folder contains scripts used throughout the project for co-folding analysis, structure processing, sequence handling, ligand sampling, and DNA optimization.

---

## Scripts overview

| Script | Purpose |
|--------|---------|
| [`check_cofolding_success.py`](#check_cofolding_successpy) | Analyze co-folding predictions to assess biosensor folding success |
| [`conda_exporter.sh`](#conda_exportersh) | Export the active Conda environment for reproducibility |
| [`extract_monomer.py`](#extract_monomerpy) | Isolate a single chain from a multimeric PDB file |
| [`extract_seq_pdb.py`](#extract_seq_pdbpy) | Extract amino acid sequence from PDB coordinates as FASTA |
| [`fasta2xlsx.py`](#fasta2xlsxpy) | Convert a FASTA file into a structured Excel table |
| [`hth_family_screening.py`](#hth_family_screeningpy) | Screen for HTH family motifs relevant to biosensor design |
| [`ligand_sampling.py`](#ligand_samplingpy) | Sample ligand poses within a protein binding pocket |
| [`protein_to_codon_opt_DNA.py`](#protein_to_codon_opt_dnapy) | Back-translate a protein sequence into codon-optimized DNA |
| [`run_boltz.py`](#run_boltzpy) | Wrapper to execute the Boltz structural prediction pipeline |

---

## Script details

### `check_cofolding_success.py`
Analyzes the output of co-folding predictions to determine whether the designed biosensor folds into its intended conformation in the presence of its target ligand.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--pred` | Path to the predicted PDB structure |
| `--ref` | Path to the reference/target PDB structure |

```bash
python check_cofolding_success.py --pred predicted.pdb --ref target.pdb
```

---

### `conda_exporter.sh`
Exports the currently active Conda environment and its dependencies to a `.yml` file, ensuring reproducibility across different machines.

> No arguments required. Operates on the active Conda environment.

```bash
bash conda_exporter.sh
```

---

### `extract_monomer.py`
Parses a multimeric or complex PDB file and isolates a single chain into a new PDB file.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input` | Path to the multi-chain PDB file |
| `--chain` | Chain identifier to extract (e.g. `A`) |
| `--output` | Path for the output monomer PDB file |

```bash
python extract_monomer.py --input complex.pdb --chain A --output monomer_A.pdb
```

---

### `extract_seq_pdb.py`
Extracts the amino acid sequence directly from the atomic coordinates of a PDB file and outputs it in FASTA format.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--pdb` | Path to the input PDB file |
| `--out` | Path for the output FASTA file |

```bash
python extract_seq_pdb.py --pdb input.pdb --out sequence.fasta
```

---

### `fasta2xlsx.py`
Converts a FASTA sequence file into a structured Excel `.xlsx` format, useful for maintaining variant libraries or inspecting sequences in a tabular format.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--fasta` | Path to the input `.fasta` or `.fa` file |
| `--excel` | Path for the output `.xlsx` file |

```bash
python fasta2xlsx.py --fasta sequences.fasta --excel sequences_output.xlsx
```

---

### `hth_family_screening.py`
Screens sequence or structural databases to identify and analyze Helix-Turn-Helix (HTH) family motifs relevant to biosensor chassis selection.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--input` | Input sequences or alignment file |
| `--output` | Path for the output CSV of HTH hits |

```bash
python hth_family_screening.py --input alignments.fasta --output hth_hits.csv
```

---

### `ligand_sampling.py`
Samples ligand poses and conformations within a defined protein binding pocket to prepare for downstream binding site design or docking evaluations.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--receptor` | Path to the receptor structure PDB file |
| `--ligand` | Path to the ligand molecule file (SDF/MOL2) |
| `--samples` | Number of poses to sample |

```bash
python ligand_sampling.py --receptor protein.pdb --ligand molecule.sdf --samples 100
```

---

### `protein_to_codon_opt_DNA.py`
Takes a protein amino acid sequence and back-translates it into a codon-optimized DNA sequence suitable for expression in a specified host organism.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--prot` | Path to the input protein FASTA file |
| `--host` | Target host organism (e.g. `e_coli`) |
| `--out` | Path for the output optimized DNA FASTA file |

```bash
python protein_to_codon_opt_DNA.py --prot target.fasta --host e_coli --out opt_dna.fasta
```

---

### `run_boltz.py`
A wrapper script to execute the Boltz structural prediction pipeline for a given set of targets or sequences.

**Arguments:**
| Argument | Description |
|----------|-------------|
| `--config` | Path to the Boltz configuration YAML file |
| `--outdir` | Directory where results will be saved |

```bash
python run_boltz.py --config boltz_config.yaml --outdir results/
```
