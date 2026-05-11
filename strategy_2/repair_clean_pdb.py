#!/usr/bin/env python3
"""
repair_clean_pdb.py

Given a structure file (PDB or CIF) this script can:
- Removes waters
- Remove hetero atoms
- Keep the chains indicated in the input
- Repair missing backbone atoms (N, C) by reconstructing them geometrically
- Detect missing residues by analyzing CA-CA distances and residue numbering

[USAGE]
python repair_clean_pdb.py \
    --input structure.pdb \
    --output repaired_structure.pdb \
    --remove_waters \
    --remove_heteroatoms \
    --keep_heteroatoms HEM,ZN \
    --keep_chains A,D \
    --repair_backbone
"""

import argparse
import numpy as np
from pathlib import Path
from collections import defaultdict
import logging

from prody import parsePDB, writePDB, AtomGroup, LOGGER
try:
    from prody.proteins import parseCIF
except ImportError:
    try:
        from prody.proteins.ciffile import parseCIF
    except ImportError:
        parseCIF = None

# Try to import biopython for CIF fallback
try:
    from Bio.PDB import MMCIFParser, PDBIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

LOGGER.verbosity = 'warning'
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def cif_to_prody_via_biopython(cif_path):
    """Convert CIF file to ProDy structure using biopython."""
    if not BIOPYTHON_AVAILABLE:
        logger.error("Biopython not available. Install with: conda install biopython")
        return None
    
    try:
        import tempfile
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', str(cif_path))
        
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            tmp_pdb = tmp.name
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_pdb)
        
        protein = parsePDB(tmp_pdb)
        Path(tmp_pdb).unlink()
        logger.info(f"Successfully converted CIF to ProDy structure using biopython")
        return protein
    except Exception as e:
        logger.error(f"Biopython CIF conversion failed: {str(e)}")
        return None


def parse_structure(structure_path):
    """Parse a protein structure file (PDB or CIF format)."""
    structure_path = Path(structure_path)
    if not structure_path.exists():
        logger.error(f"File not found: {structure_path}")
        return None
    
    suffix = structure_path.suffix.lower()
    try:
        if suffix == '.cif':
            logger.info(f"Parsing CIF file: {structure_path}")
            if parseCIF is not None:
                return parseCIF(str(structure_path))
            logger.info("ProDy CIF parser not available, attempting biopython conversion...")
            return cif_to_prody_via_biopython(str(structure_path))
        elif suffix == '.pdb':
            logger.info(f"Parsing PDB file: {structure_path}")
            return parsePDB(str(structure_path))
        else:
            logger.error(f"Unsupported file format: {suffix}. Use .pdb or .cif")
            return None
    except Exception as e:
        logger.error(f"Error parsing {structure_path}: {str(e)}")
        return None


def clean_structure(protein, remove_waters, remove_heteroatoms, keep_heteroatoms, keep_chains):
    """Apply ProDy selections to clean the structure."""
    selections = []

    if keep_chains:
        chain_list = keep_chains.replace(',', ' ')
        selections.append(f"(chain {chain_list})")
        
    if remove_waters:
        selections.append("(not water)")
        
    if remove_heteroatoms:
        if keep_heteroatoms:
            het_list = keep_heteroatoms.replace(',', ' ')
            selections.append(f"(not hetero or resname {het_list})")
        else:
            selections.append("(not hetero)")

    if selections:
        sel_string = " and ".join(selections)
        logger.info(f"Applying cleaning selection: {sel_string}")
        protein = protein.select(sel_string)
        if protein is not None:
            protein = protein.copy()
        else:
            logger.error("Cleaning resulted in an empty structure! Check your flags.")
    
    return protein


def repair_missing_backbone_atoms(protein):
    """Repair missing backbone atoms (N or C) by reconstructing them geometrically."""
    chains = set(protein.getChids())
    atoms_to_add = []
    current_max_serial = max(protein.getSerials()) if len(protein.getSerials()) > 0 else 0
    
    for chain_id in chains:
        chain_resnums = protein.select(f'chain {chain_id}').getResnums()
        chain_resnums = sorted(set(chain_resnums))
        
        for resnum in chain_resnums:
            residue = protein.select(f'chain {chain_id} and resnum {resnum}')
            if residue is None or len(residue) == 0:
                continue
            
            resname = residue.getResnames()[0]
            atom_names = set(residue.getNames())
            
            # Check if N is missing
            if 'N' not in atom_names and 'CA' in atom_names:
                try:
                    ca_coords = residue.select('name CA').getCoords()[0]
                    prev_resnum = resnum - 1
                    prev_residue = protein.select(f'chain {chain_id} and resnum {prev_resnum} and name C')
                    if prev_residue is not None and len(prev_residue) > 0:
                        prev_c = prev_residue.getCoords()[0]
                        direction = (ca_coords - prev_c) / np.linalg.norm(ca_coords - prev_c)
                        n_coord = ca_coords - 1.45 * direction
                    else:
                        n_coord = ca_coords + np.array([-0.5, 0.0, 0.0])
                    
                    atoms_to_add.append({
                        'name': 'N', 'resnum': resnum, 'resname': resname,
                        'chain': chain_id, 'coords': n_coord, 'bfactor': 0.0
                    })
                    logger.info(f"Reconstructed missing N atom for {chain_id}{resnum}")
                except Exception as e:
                    logger.warning(f"Could not reconstruct N for {chain_id}{resnum}: {str(e)}")
            
            # Check if C is missing
            if 'C' not in atom_names and 'CA' in atom_names:
                try:
                    ca_coords = residue.select('name CA').getCoords()[0]
                    next_resnum = resnum + 1
                    next_residue = protein.select(f'chain {chain_id} and resnum {next_resnum} and name N')
                    if next_residue is not None and len(next_residue) > 0:
                        next_n = next_residue.getCoords()[0]
                        direction = (next_n - ca_coords) / np.linalg.norm(next_n - ca_coords)
                        c_coord = ca_coords + 1.52 * direction
                    else:
                        c_coord = ca_coords + np.array([0.5, 0.0, 0.0])
                    
                    atoms_to_add.append({
                        'name': 'C', 'resnum': resnum, 'resname': resname,
                        'chain': chain_id, 'coords': c_coord, 'bfactor': 0.0
                    })
                    logger.info(f"Reconstructed missing C atom for {chain_id}{resnum}")
                except Exception as e:
                    logger.warning(f"Could not reconstruct C for {chain_id}{resnum}: {str(e)}")
    
    if atoms_to_add:
        logger.info(f"Adding {len(atoms_to_add)} reconstructed backbone atoms")
        new_atoms = AtomGroup('Reconstructed Backbone')
        
        new_atoms.setSerials(np.array([current_max_serial + i + 1 for i in range(len(atoms_to_add))]))
        new_atoms.setNames(np.array([a['name'] for a in atoms_to_add]))
        new_atoms.setResnames(np.array([a['resname'] for a in atoms_to_add]))
        new_atoms.setResnums(np.array([a['resnum'] for a in atoms_to_add]))
        new_atoms.setChids(np.array([a['chain'] for a in atoms_to_add]))
        new_atoms.setCoords(np.array([a['coords'] for a in atoms_to_add]))
        new_atoms.setBetas(np.array([a['bfactor'] for a in atoms_to_add]))
        
        protein = protein + new_atoms
    
    return protein


def detect_chain_breaks(protein, gap_threshold_distance=4.0, gap_threshold_resnum=1):
    """Detect missing residues by analyzing CA-CA distances and residue numbering."""
    gaps = defaultdict(list)
    ca_atoms = protein.select('name CA')
    if ca_atoms is None or len(ca_atoms) == 0:
        logger.warning("No CA atoms found in protein, skipping gap detection")
        return gaps
    
    chains = defaultdict(list)
    for atom in ca_atoms:
        chains[atom.getChid()].append(atom)
    
    for chain_id, atoms in chains.items():
        atoms_sorted = sorted(atoms, key=lambda a: a.getResnum())
        for i in range(len(atoms_sorted) - 1):
            atom1, atom2 = atoms_sorted[i], atoms_sorted[i + 1]
            resnum1, resnum2 = atom1.getResnum(), atom2.getResnum()
            
            resnum_gap = resnum2 - resnum1
            distance = np.linalg.norm(atom2.getCoords() - atom1.getCoords())
            
            if resnum_gap > gap_threshold_resnum or distance > gap_threshold_distance:
                missing_start, missing_end = resnum1 + 1, resnum2 - 1
                if missing_start <= missing_end:
                    gaps[chain_id].append((missing_start, missing_end, resnum1, resnum2))
                    logger.info(
                        f"Chain {chain_id}: Gap detected between residues {resnum1} and {resnum2} "
                        f"(distance: {distance:.2f}Å, resnum_gap: {resnum_gap})"
                    )
    return dict(gaps)


def fix_occupancy_in_pdb(pdb_path):
    """Fix zero/partial occupancy values in PDB file output by biopython."""
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    fixed_lines = []
    occupancy_fixed = 0
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            try:
                occupancy = float(line[54:60])
                if occupancy < 1.0:
                    fixed_occupancy = " " * 2 + "1.00"
                    line = line[:54] + fixed_occupancy + line[60:]
                    occupancy_fixed += 1
            except (ValueError, IndexError):
                pass
        fixed_lines.append(line)
    
    if occupancy_fixed > 0:
        with open(pdb_path, 'w') as f:
            f.writelines(fixed_lines)
        logger.info(f"Fixed occupancy for {occupancy_fixed} atoms.")


def main():
    parser = argparse.ArgumentParser(description="Clean and repair PDB/CIF structure files.")
    parser.add_argument('--input', required=True, help='Path to input structure file (.pdb or .cif)')
    parser.add_argument('--output', required=True, help='Path for output repaired structure (.pdb)')
    parser.add_argument('--remove_waters', action='store_true', help='Remove water molecules')
    parser.add_argument('--remove_heteroatoms', action='store_true', help='Remove heteroatoms (HETATM)')
    parser.add_argument('--keep_heteroatoms', type=str, help='Comma-separated list of heteroatoms to keep (e.g., HEM,ZN)')
    parser.add_argument('--keep_chains', type=str, help='Comma-separated list of chains to keep (e.g., A,D)')
    parser.add_argument('--repair_backbone', action='store_true', help='Repair missing N and C backbone atoms')
    
    args = parser.parse_args()
    
    # 1. Parse Structure
    protein = parse_structure(args.input)
    if protein is None:
        logger.error("Failed to parse structure. Exiting.")
        return
        
    # 2. Clean Structure (waters, heteroatoms, chains)
    if any([args.remove_waters, args.remove_heteroatoms, args.keep_chains]):
        protein = clean_structure(
            protein, 
            remove_waters=args.remove_waters, 
            remove_heteroatoms=args.remove_heteroatoms, 
            keep_heteroatoms=args.keep_heteroatoms, 
            keep_chains=args.keep_chains
        )
        if protein is None:
            return

    # 3. Backbone Repair
    if args.repair_backbone:
        logger.info("Initiating backbone repair...")
        protein = repair_missing_backbone_atoms(protein)

    # 4. Chain Break Detection
    logger.info("Running chain break detection...")
    gaps = detect_chain_breaks(protein)
    if not gaps:
        logger.info("No chain breaks detected.")
        
    # 5. Write Output
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    writePDB(str(output_path), protein)
    logger.info(f"Successfully wrote output to {output_path}")
    
    # Optional: Fix occupancies if CIF > PDB mangled them
    fix_occupancy_in_pdb(str(output_path))


if __name__ == "__main__":
    main()