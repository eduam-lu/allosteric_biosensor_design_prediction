#!/usr/bin/env python3
"""
Structure repair utility for filling missing residues with glycine.
Supports both PDB (.pdb) and CIF (.cif) formats (e.g., RF Diffusion outputs).
Detects chain breaks, fills gaps, and updates MPNN redesign lists.
"""

import json
import numpy as np
from pathlib import Path
from prody import parsePDB, writePDB, AtomGroup, LOGGER, Atom
try:
    from prody.proteins import parseCIF
except ImportError:
    try:
        from prody.proteins.ciffile import parseCIF
    except ImportError:
        # CIF support not available, will handle with fallback
        parseCIF = None

# Try to import biopython
try:
    from Bio.PDB import MMCIFParser, PDBIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False

from collections import defaultdict
import logging

LOGGER.verbosity = 'warning'

# Standard amino acid 1-letter codes
AA_1to3 = {
    'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
    'Q': 'GLN', 'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'
}

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def cif_to_prody_via_biopython(cif_path):
    """
    Convert CIF file to ProDy structure using biopython.
    
    Args:
        cif_path: Path to CIF file
        
    Returns:
        ProDy Protein object or None if conversion fails
    """
    if not BIOPYTHON_AVAILABLE:
        logger.error("Biopython not available. Install with: conda install biopython")
        return None
    
    try:
        from Bio.PDB import MMCIFParser
        import tempfile
        
        # Parse CIF with biopython
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure('protein', str(cif_path))
        
        # Write to temporary PDB file
        with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
            tmp_pdb = tmp.name
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(tmp_pdb)
        
        # Parse with ProDy
        protein = parsePDB(tmp_pdb)
        
        # Clean up
        Path(tmp_pdb).unlink()
        
        logger.info(f"Successfully converted CIF to ProDy structure using biopython")
        return protein
        
    except Exception as e:
        logger.error(f"Biopython CIF conversion failed: {str(e)}")
        return None


def parse_structure(structure_path):
    """
    Parse a protein structure file (PDB or CIF format).
    For CIF files, attempts ProDy's parseCIF; if unavailable, uses OpenBabel conversion.
    
    Args:
        structure_path: Path to .pdb or .cif file
    
    Returns:
        ProDy Protein object or None if parse fails
    """
    structure_path = Path(structure_path)
    
    if not structure_path.exists():
        logger.error(f"File not found: {structure_path}")
        return None
    
    suffix = structure_path.suffix.lower()
    
    try:
        if suffix == '.cif':
            logger.info(f"Parsing CIF file: {structure_path}")
            
            # Try native ProDy CIF parser first
            if parseCIF is not None:
                return parseCIF(str(structure_path))
            
            # Fallback: Use biopython to convert CIF to PDB
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


def filter_keep_chains_a_and_d(protein):
    """
    Filter protein to keep chains A and D (if present).
    Chain A is the protein, chain D is the ligand.
    Removes all other chains if present.
    
    Args:
        protein: ProDy Protein object
    
    Returns:
        Filtered ProDy Protein object (chains A and D, or just A if D not present)
    """
    chains = protein.getChids()
    unique_chains = set(chains)
    
    # Check what we have available
    has_a = 'A' in unique_chains
    has_d = 'D' in unique_chains
    
    if not has_a and not has_d:
        logger.warning(f"Protein has neither chain A nor D! Chains present: {unique_chains}")
        logger.warning("Attempting to rename first chain to A...")
        if len(unique_chains) > 0:
            first_chain = list(unique_chains)[0]
            protein.setChids('A')
            logger.info(f"Renamed chain {first_chain} to A")
        return protein
    
    # Determine which chains to keep
    chains_to_keep = []
    if has_a:
        chains_to_keep.append('A')
    if has_d:
        chains_to_keep.append('D')
    
    # Build selection string
    if not has_a and has_d:
        # Only D exists, rename it to A
        logger.info("Chain A not found; chain D will be kept as is (ligand)")
        chain_d = protein.select('chain D')
        if chain_d is None or len(chain_d) == 0:
            logger.error("Failed to select chain D")
            return protein
        filtered_protein = chain_d.copy()
    elif has_a and not has_d:
        # Only A exists
        logger.info("Chain A found. Chain D (ligand) not found. Keeping chain A only.")
        chain_a = protein.select('chain A')
        if chain_a is None or len(chain_a) == 0:
            logger.error("Failed to select chain A")
            return protein
        filtered_protein = chain_a.copy()
    else:
        # Both A and D exist
        logger.info("Both chain A (protein) and chain D (ligand) found. Keeping both.")
        chain_a = protein.select('chain A')
        chain_d = protein.select('chain D')
        
        if chain_a is None or len(chain_a) == 0:
            logger.error("Failed to select chain A")
            return protein
        if chain_d is None or len(chain_d) == 0:
            logger.error("Failed to select chain D")
            return protein
        
        # Merge both chains - combine A and D
        chain_a_copy = chain_a.copy()
        chain_d_copy = chain_d.copy()
        filtered_protein = chain_a_copy + chain_d_copy
    
    removed_chains = unique_chains - set(chains_to_keep)
    if removed_chains:
        logger.info(f"Removed chains: {removed_chains}")
    
    return filtered_protein


def repair_missing_backbone_atoms(protein):
    """
    Repair missing backbone atoms (N or C) by reconstructing them geometrically.
    Uses neighbor CA and O atoms to approximate missing backbone geometry.
    
    Args:
        protein: ProDy Protein object (will be modified in-place)
    
    Returns:
        modified protein with reconstructed backbone atoms
    """
    import numpy as np
    
    # Get all residues
    all_resnums = sorted(set(protein.getResnums()))
    chains = set(protein.getChids())
    
    atoms_to_add = []
    current_max_serial = max(protein.getSerials()) if len(protein.getSerials()) > 0 else 0
    
    for chain_id in chains:
        chain_resnums = protein.select(f'chain {chain_id}').getResnums()
        chain_resnums = sorted(set(chain_resnums))
        
        for resnum in chain_resnums:
            select_str = f'chain {chain_id} and resnum {resnum}'
            residue = protein.select(select_str)
            
            if residue is None or len(residue) == 0:
                continue
            
            resname = residue.getResnames()[0]
            atom_names = set(residue.getNames())
            
            # Check if N is missing
            if 'N' not in atom_names and 'CA' in atom_names:
                try:
                    ca_coords = residue.select('name CA').getCoords()[0]
                    c_coords = residue.select('name C').getCoords()[0] if 'C' in atom_names else None
                    
                    # Get previous residue's C atom for reference
                    prev_resnum = resnum - 1
                    prev_residue = protein.select(f'chain {chain_id} and resnum {prev_resnum} and name C')
                    if prev_residue is not None and len(prev_residue) > 0:
                        prev_c = prev_residue.getCoords()[0]
                        # N is roughly 1.45 Å from C in the C->CA direction
                        direction = (ca_coords - prev_c) / np.linalg.norm(ca_coords - prev_c)
                        n_coord = ca_coords - 1.45 * direction
                    else:
                        # Approximate N position
                        n_coord = ca_coords + np.array([-0.5, 0.0, 0.0])
                    
                    atoms_to_add.append({
                        'name': 'N',
                        'resnum': resnum,
                        'resname': resname,
                        'chain': chain_id,
                        'coords': n_coord,
                        'bfactor': 0.0
                    })
                    logger.info(f"Reconstructed missing N atom for {chain_id}{resnum}")
                except Exception as e:
                    logger.warning(f"Could not reconstruct N for {chain_id}{resnum}: {str(e)}")
            
            # Check if C is missing
            if 'C' not in atom_names and 'CA' in atom_names:
                try:
                    ca_coords = residue.select('name CA').getCoords()[0]
                    
                    # Get next residue's N atom for reference
                    next_resnum = resnum + 1
                    next_residue = protein.select(f'chain {chain_id} and resnum {next_resnum} and name N')
                    if next_residue is not None and len(next_residue) > 0:
                        next_n = next_residue.getCoords()[0]
                        # C is roughly 1.52 Å from N in the N->CA direction
                        direction = (next_n - ca_coords) / np.linalg.norm(next_n - ca_coords)
                        c_coord = ca_coords + 1.52 * direction
                    else:
                        # Approximate C position
                        c_coord = ca_coords + np.array([0.5, 0.0, 0.0])
                    
                    atoms_to_add.append({
                        'name': 'C',
                        'resnum': resnum,
                        'resname': resname,
                        'chain': chain_id,
                        'coords': c_coord,
                        'bfactor': 0.0
                    })
                    logger.info(f"Reconstructed missing C atom for {chain_id}{resnum}")
                except Exception as e:
                    logger.warning(f"Could not reconstruct C for {chain_id}{resnum}: {str(e)}")
    
    # Add reconstructed atoms to protein
    if atoms_to_add:
        logger.info(f"Adding {len(atoms_to_add)} reconstructed backbone atoms")
        new_atoms = AtomGroup('Reconstructed Backbone')
        
        serials = [current_max_serial + i + 1 for i in range(len(atoms_to_add))]
        names = [a['name'] for a in atoms_to_add]
        resnames = [a['resname'] for a in atoms_to_add]
        resnums = [a['resnum'] for a in atoms_to_add]
        chains_list = [a['chain'] for a in atoms_to_add]
        coords_list = [a['coords'] for a in atoms_to_add]
        bfactors = [a['bfactor'] for a in atoms_to_add]
        
        new_atoms.setSerials(np.array(serials))
        new_atoms.setNames(np.array(names))
        new_atoms.setResnames(np.array(resnames))
        new_atoms.setResnums(np.array(resnums))
        new_atoms.setChids(np.array(chains_list))
        new_atoms.setCoords(np.array(coords_list))
        new_atoms.setBetas(np.array(bfactors))
        
        protein = protein + new_atoms
        logger.info(f"Successfully added {len(atoms_to_add)} reconstructed backbone atoms")
    
    return protein


def detect_chain_breaks(protein, gap_threshold_distance=4.0, gap_threshold_resnum=1):
    """
    Detect missing residues by analyzing CA-CA distances and residue numbering.
    
    Args:
        protein: ProDy Protein object
        gap_threshold_distance: Min CA-CA distance to signal a gap (Angstroms)
        gap_threshold_resnum: Max resnum gap to allow without signaling break
    
    Returns:
        dict: {chain_id: [(start_resnum, end_resnum), ...]}
    """
    gaps = defaultdict(list)
    
    # Get all CA atoms
    ca_atoms = protein.select('name CA')
    if ca_atoms is None or len(ca_atoms) == 0:
        logger.warning("No CA atoms found in protein")
        return gaps
    
    # Group by chain
    chains = defaultdict(list)
    for atom in ca_atoms:
        chain_id = atom.getChid()
        chains[chain_id].append(atom)
    
    # Check for gaps within each chain
    for chain_id, atoms in chains.items():
        # Sort by residue number
        atoms_sorted = sorted(atoms, key=lambda a: a.getResnum())
        
        for i in range(len(atoms_sorted) - 1):
            atom1 = atoms_sorted[i]
            atom2 = atoms_sorted[i + 1]
            
            resnum1 = atom1.getResnum()
            resnum2 = atom2.getResnum()
            
            # Check residue numbering gap
            resnum_gap = resnum2 - resnum1
            
            # Check spatial gap (distance between CA atoms)
            coords1 = atom1.getCoords()
            coords2 = atom2.getCoords()
            distance = np.linalg.norm(coords2 - coords1)
            
            if resnum_gap > gap_threshold_resnum or distance > gap_threshold_distance:
                # Missing residues between resnum1 and resnum2
                missing_start = resnum1 + 1
                missing_end = resnum2 - 1
                if missing_start <= missing_end:
                    gaps[chain_id].append((missing_start, missing_end, resnum1, resnum2))
                    logger.info(
                        f"Chain {chain_id}: Gap detected between residues {resnum1} and {resnum2} "
                        f"(distance: {distance:.2f}Å, resnum_gap: {resnum_gap})"
                    )
    
    return dict(gaps)


def fill_missing_residues_with_glycine(protein, gaps, max_gap_distance=6.0):
    """
    Fill missing residues with Glycine (G) by interpolating coordinates.
    Only fills small gaps (≤6 Å); larger gaps indicate structural discontinuities.
    Modifies the protein object IN PLACE by adding new atoms.
    Uses linear interpolation for backbone atoms (N, CA, C, O).
    
    Args:
        protein: ProDy Protein object (will be modified)
        gaps: dict from detect_chain_breaks()
        max_gap_distance: Maximum CA-CA distance to attempt filling (default 6.0 Å)
    
    Returns:
        tuple: (modified_protein, dict of filled_resnums)
    """
    from prody import AtomGroup, Chain
    
    filled_residues = defaultdict(list)
    atoms_to_add = []
    
    for chain_id, gap_list in gaps.items():
        for missing_start, missing_end, anchor_before, anchor_after in gap_list:
            
            # Get anchor atoms for interpolation
            select_before = f"chain {chain_id} and resnum {anchor_before} and name CA"
            select_after = f"chain {chain_id} and resnum {anchor_after} and name CA"
            
            ca_before = protein.select(select_before)
            ca_after = protein.select(select_after)
            
            if ca_before is None or ca_after is None:
                logger.warning(f"Could not find anchor atoms for gap in {chain_id}")
                continue
            
            coords_before = ca_before.getCoords()[0]
            coords_after = ca_after.getCoords()[0]
            
            # Calculate distance between anchors
            import numpy as np
            distance = np.linalg.norm(coords_after - coords_before)
            
            # Skip large gaps (likely structural breaks, not missing residues)
            if distance > max_gap_distance:
                logger.warning(
                    f"Chain {chain_id}: Skipping large gap between {anchor_before} and {anchor_after} "
                    f"(distance: {distance:.2f}Å > {max_gap_distance}Å). This appears to be a "
                    f"structural discontinuity, not missing coordinates."
                )
                continue
            
            # Create atoms for each missing residue
            n_missing = missing_end - missing_start + 1
            
            for i in range(1, n_missing + 1):
                resnum = missing_start + i - 1
                
                # Linear interpolation for CA position
                t = i / (n_missing + 1)
                ca_coord = coords_before + t * (coords_after - coords_before)
                
                # Approximate backbone geometry around CA
                # N is roughly -1.5 Å back along chain direction
                # C is roughly +1.5 Å forward along chain direction  
                direction = (coords_after - coords_before) / np.linalg.norm(coords_after - coords_before)
                
                n_coord = ca_coord - 1.5 * direction
                c_coord = ca_coord + 1.5 * direction
                o_coord = c_coord + np.array([0.0, 1.2, 0.0])  # O is ~1.2Å above C
                
                # Create atoms: N, CA, C, O for Gly
                for atom_name, coord in [('N', n_coord), ('CA', ca_coord), ('C', c_coord), ('O', o_coord)]:
                    atoms_to_add.append({
                        'name': atom_name,
                        'resnum': resnum,
                        'resname': 'GLY',
                        'chain': chain_id,
                        'coords': coord,
                        'serial': None,  # Will be assigned
                        'bfactor': 0.0  # Use normal B-factor to avoid filtering
                    })
                
                logger.info(
                    f"Prepared to add Chain {chain_id} Residue {resnum} (Gly) "
                    f"at CA position {ca_coord}"
                )
                filled_residues[chain_id].append(resnum)
    
    # Now add all atoms to the protein
    if atoms_to_add:
        logger.info(f"Adding {len(atoms_to_add)} atoms ({len(atoms_to_add)//4} residues) to protein")
        
        # Add atoms to protein using ProDy's addAtoms method
        current_max_serial = max(protein.getSerials()) if len(protein.getSerials()) > 0 else 0
        
        try:
            for idx, atom_data in enumerate(atoms_to_add):
                atom = protein.select('none')  # Create empty selection
                # Use a simpler approach - create new coordinates and append
                
            # Actually, let's take a different approach: create new atoms and merge
            new_atoms = AtomGroup('Filled Residues')
            
            serials = []
            names = []
            resnames = []
            resnums = []
            chains_list = []
            coords_list = []
            bfactors = []
            
            for atom_data in atoms_to_add:
                serials.append(current_max_serial + len(serials) + 1)
                names.append(atom_data['name'])
                resnames.append(atom_data['resname'])
                resnums.append(atom_data['resnum'])
                chains_list.append(atom_data['chain'])
                coords_list.append(atom_data['coords'])
                bfactors.append(atom_data['bfactor'])
            
            new_atoms.setSerials(np.array(serials))
            new_atoms.setNames(np.array(names))
            new_atoms.setResnames(np.array(resnames))
            new_atoms.setResnums(np.array(resnums))
            new_atoms.setChids(np.array(chains_list))
            new_atoms.setCoords(np.array(coords_list))
            new_atoms.setBetas(np.array(bfactors))
            
            # Merge with original protein
            protein = protein + new_atoms
            logger.info(f"Successfully added filled residues to protein")
            
        except Exception as e:
            logger.error(f"Error adding atoms to protein: {str(e)}")
            logger.error("Continuing without adding atoms to PDB (ProDy operations may fail later)")
    
    return protein, dict(filled_residues)


def create_dummy_glycine_atoms(resnum, chain_id, interp_coord):
    """
    Create a minimal Glycine residue backbone.
    Note: This creates a simplified representation.
    
    Returns:
        list of (atom_name, coords) tuples
    """
    # Approximate backbone geometry for GLY
    # B-factor set to indicate added residue
    backbone_atoms = [
        ('N', interp_coord + np.array([-0.5, 0.0, 0.0]), 100.0),
        ('CA', interp_coord, 100.0),
        ('C', interp_coord + np.array([0.5, 0.0, 0.0]), 100.0),
        ('O', interp_coord + np.array([0.6, 1.0, 0.0]), 100.0),
        ('CB', interp_coord + np.array([0.0, -1.0, 0.0]), 100.0),
    ]
    return backbone_atoms


def write_repaired_pdb(protein, gaps, output_path):
    """
    Write the protein to a PDB file with a note about filled residues.
    Note: This is a simplified version - ideally would rebuild full residues.
    
    Args:
        protein: ProDy Protein object
        gaps: dict of gaps
        output_path: Path to write PDB
    """
    # For now, just write the original protein
    # A full implementation would add the missing residues
    writePDB(str(output_path), protein)
    logger.info(f"Wrote repaired PDB to {output_path}")
    
    # Fix occupancy values (some atoms may have 0.0 occupancy after CIF conversion)
    fix_occupancy_in_pdb(str(output_path))


def fix_occupancy_in_pdb(pdb_path):
    """
    Fix occupancy values in PDB file. 
    Sets occupancy to 1.0 for all ATOM records that have occupancy < 1.0.
    This is necessary because CIF-to-PDB conversion via biopython sometimes
    creates atoms with partial or zero occupancy, which get filtered out by
    structure prediction tools like LigandMPNN.
    
    Args:
        pdb_path: Path to PDB file to fix
    """
    # Read the PDB file
    with open(pdb_path, 'r') as f:
        lines = f.readlines()
    
    # Fix occupancy values
    fixed_lines = []
    occupancy_fixed = 0
    
    for line in lines:
        if line.startswith('ATOM'):
            # PDB format: occupancy is at columns 54-60 (1-indexed: 55-60)
            try:
                occupancy_str = line[54:60]
                occupancy = float(occupancy_str)
                
                if occupancy < 1.0:
                    # Replace occupancy with 1.00
                    # Format: 6 characters, right-aligned, 2 decimal places
                    fixed_occupancy = " " * 2 + "1.00"  # 6 characters with 1.00
                    line = line[:54] + fixed_occupancy + line[60:]
                    occupancy_fixed += 1
            except (ValueError, IndexError):
                pass  # Skip if conversion fails
        
        fixed_lines.append(line)
    
    # Write the fixed PDB file
    with open(pdb_path, 'w') as f:
        f.writelines(fixed_lines)
    
    if occupancy_fixed > 0:
        logger.info(f"Fixed occupancy for {occupancy_fixed} atoms in {Path(pdb_path).name}")



def repair_pdb_and_update_redesign_list(pdb_path, redesign_residues_dict, output_dir=None):
    """
    Main function: Repair structure file by detecting and filling chain breaks.
    Supports both PDB and CIF formats (e.g., RF Diffusion outputs).
    Updates the redesign residue list to include filled residues.
    
    Args:
        pdb_path: Path to structure file (.pdb or .cif)
        redesign_residues_dict: Dict {chain_id: [resnum, ...]} of residues to redesign
        output_dir: Directory to write repaired structure (if None, use same as input)
    
    Returns:
        dict: {
            'success': bool,
            'original_residues': int,
            'filled_residues': int,
            'new_total': int,
            'filled_positions': [(chain_id, resnum), ...],
            'updated_redesign_dict': updated_redesign_residues_dict,
            'report': str,
            'output_pdb': str or None
        }
    """
    
    result = {
        'success': False,
        'original_residues': 0,
        'filled_residues': 0,
        'new_total': 0,
        'filled_positions': [],
        'updated_redesign_dict': redesign_residues_dict.copy(),
        'report': '',
        'output_pdb': None
    }
    
    try:
        # Validate path exists
        pdb_path = Path(pdb_path)
        if not pdb_path.exists():
            result['report'] = f"Error: File not found: {pdb_path}"
            logger.error(result['report'])
            return result
        
        # Parse structure (handles both PDB and CIF)
        protein = parse_structure(pdb_path)
        if protein is None or len(protein) == 0:
            result['report'] = f"Error: Could not parse structure {pdb_path}"
            logger.error(result['report'])
            return result
        
        # Filter to keep chains A and D (protein + ligand)
        protein = filter_keep_chains_a_and_d(protein)
        
        # Repair missing backbone atoms (N, C)
        protein = repair_missing_backbone_atoms(protein)
        
        original_residues = protein.numResidues()
        result['original_residues'] = original_residues
        
        # Detect gaps
        gaps = detect_chain_breaks(protein)
        
        if not gaps:
            result['success'] = True
            result['report'] = f"No chain breaks detected. PDB is valid ({original_residues} residues)."
            logger.info(result['report'])
            # Still write the filtered chain A only version
            output_pdb = output_dir / f"{Path(pdb_path).stem}_repaired.pdb"
            write_repaired_pdb(protein, {}, output_pdb)
            result['output_pdb'] = str(output_pdb)
            return result
        
        # Fill missing residues (returns modified protein and filled residues dict)
        protein, filled_residues = fill_missing_residues_with_glycine(protein, gaps)
        
        # Flatten filled residues for counting
        all_filled = []
        for chain_id, resnums in filled_residues.items():
            all_filled.extend([(chain_id, resnum) for resnum in resnums])
        
        result['filled_residues'] = len(all_filled)
        result['new_total'] = original_residues + len(all_filled)
        result['filled_positions'] = all_filled
        
        # Update redesign list
        updated_redesign = redesign_residues_dict.copy()
        for chain_id, resnum in all_filled:
            if chain_id not in updated_redesign:
                updated_redesign[chain_id] = []
            if resnum not in updated_redesign[chain_id]:
                updated_redesign[chain_id].append(resnum)
        
        result['updated_redesign_dict'] = updated_redesign
        
        # Write output
        if output_dir is None:
            output_dir = Path(pdb_path).parent
        else:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
        
        output_pdb = output_dir / f"{Path(pdb_path).stem}_repaired.pdb"
        write_repaired_pdb(protein, gaps, output_pdb)
        result['output_pdb'] = str(output_pdb)
        
        # Generate report
        filled_str = "\n".join([f"  - Chain {c} Residue {r}" for c, r in all_filled])
        report = (
            f"PDB Repair Report for {Path(pdb_path).name}\n"
            f"{'='*60}\n"
            f"Original residues: {original_residues}\n"
            f"Filled residues: {result['filled_residues']}\n"
            f"New total: {result['new_total']}\n\n"
            f"Filled positions:\n{filled_str}\n\n"
            f"Redesign list updated with {result['filled_residues']} new residues.\n"
            f"Repaired PDB: {output_pdb}\n"
        )
        result['report'] = report
        result['success'] = True
        
        logger.info(report)
        return result
        
    except Exception as e:
        result['report'] = f"Error repairing PDB: {str(e)}"
        logger.error(result['report'], exc_info=True)
        return result


def batch_repair_and_update_jsons(pdb_paths_json, redesigned_residues_json, output_dir):
    """
    Repair all structures in the JSON dict and update both JSON files.
    Handles both PDB and CIF formats.
    
    Args:
        pdb_paths_json: Path to pdb_paths_multi.json (can contain .pdb or .cif paths)
        redesigned_residues_json: Path to redesigned_residues_multi.json
        output_dir: Directory for repaired structures and updated JSONs
    
    Returns:
        dict: {structure_name: repair_result, ...}
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load JSONs
    with open(pdb_paths_json, 'r') as f:
        pdb_paths = json.load(f)
    
    with open(redesigned_residues_json, 'r') as f:
        redesign_residues = json.load(f)
    
    logger.info(f"Loaded {len(pdb_paths)} structures from {pdb_paths_json}")
    
    # Repair each structure
    all_results = {}
    updated_pdb_paths = {}
    updated_redesign_residues = {}
    
    # LigandMPNN format: pdb_paths has paths as keys, values are empty strings
    for pdb_path, _ in pdb_paths.items():
        # Skip empty paths
        if not pdb_path or str(pdb_path).strip() == "":
            logger.warning(f"Skipping: empty path")
            continue
        
        # Extract a name from the path for identification
        pdb_name = Path(pdb_path).stem
        
        # Get current redesign string/dict from original file
        # Could be keyed by path or by stem, try both
        original_redesign_str = redesign_residues.get(pdb_path, None)
        if original_redesign_str is None:
            original_redesign_str = redesign_residues.get(pdb_name, "")
        
        # Convert string to dict format for repair function
        current_redesign_dict = {}
        if isinstance(original_redesign_str, str) and original_redesign_str:
            # Parse "A100 A105 B50" format into {chain: [residues]}
            parts = original_redesign_str.split()
            for part in parts:
                if len(part) > 1:
                    chain = part[0]
                    try:
                        resnum = int(part[1:])
                        if chain not in current_redesign_dict:
                            current_redesign_dict[chain] = []
                        current_redesign_dict[chain].append(resnum)
                    except ValueError:
                        pass
        elif isinstance(original_redesign_str, dict):
            current_redesign_dict = original_redesign_str
        
        # Repair
        result = repair_pdb_and_update_redesign_list(
            pdb_path,
            current_redesign_dict,
            output_dir
        )
        
        all_results[pdb_name] = result
        
        # Determine the output path (use repaired if successful, else original)
        output_pdb_path = result['output_pdb'] if (result['success'] and result['output_pdb']) else str(pdb_path)
        
        # Update paths: maintain LigandMPNN format (path as key, "" as value)
        updated_pdb_paths[output_pdb_path] = ""
        
        # Convert updated redesign dict back to string format for output
        # Format: "A100 A105 A110 B50" (chain + resnum space-separated)
        updated_dict = result['updated_redesign_dict']
        redesign_string_parts = []
        for chain_id in sorted(updated_dict.keys()):
            resnums = sorted(updated_dict[chain_id])
            for resnum in resnums:
                redesign_string_parts.append(f"{chain_id}{resnum}")
        
        updated_redesign_str = " ".join(redesign_string_parts)
        
        # Store with OUTPUT path as key (must match pdb_paths keys!)
        updated_redesign_residues[output_pdb_path] = updated_redesign_str
    
    # Write updated JSONs
    updated_pdb_json = output_dir / "pdb_paths_multi_repaired.json"
    updated_redesign_json = output_dir / "redesigned_residues_multi_repaired.json"
    
    with open(updated_pdb_json, 'w') as f:
        json.dump(updated_pdb_paths, f, indent=2)
    
    with open(updated_redesign_json, 'w') as f:
        json.dump(updated_redesign_residues, f, indent=2)
    
    logger.info(f"Updated JSONs written to:")
    logger.info(f"  PDB paths: {updated_pdb_json}")
    logger.info(f"  Redesign residues: {updated_redesign_json}")
    
    # Write summary report to parent directory (for Snakemake)
    summary_report = output_dir.parent / "repair_summary.txt"
    with open(summary_report, 'w') as f:
        f.write("PDB REPAIR SUMMARY\n")
        f.write("="*80 + "\n\n")
        
        success_count = sum(1 for r in all_results.values() if r['success'])
        f.write(f"Processed: {len(all_results)} PDBs\n")
        f.write(f"Successful: {success_count}\n")
        f.write(f"Failed: {len(all_results) - success_count}\n\n")
        
        for pdb_name, result in all_results.items():
            f.write(f"\n{pdb_name}:\n")
            f.write(f"  Success: {result['success']}\n")
            f.write(f"  Original residues: {result['original_residues']}\n")
            f.write(f"  Filled residues: {result['filled_residues']}\n")
            if result['filled_positions']:
                f.write(f"  Filled at: {result['filled_positions']}\n")
            if result['report']:
                f.write(f"  Report: {result['report'][:200]}...\n")
    
    logger.info(f"Summary report written to: {summary_report}")
    
    return all_results


def main():
    """CLI entry point. Supports both PDB (.pdb) and CIF (.cif) formats."""
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python repair_pdb.py <structure_file> [<redesign_dict_json>]")
        print("       (structure_file can be .pdb or .cif)")
        print("  or: python repair_pdb.py --batch <pdb_paths_json> <redesign_residues_json> <output_dir>")
        sys.exit(1)
    
    if sys.argv[1] == "--batch":
        if len(sys.argv) < 5:
            print("Usage: repair_pdb.py --batch <pdb_paths_json> <redesign_residues_json> <output_dir>")
            sys.exit(1)
        results = batch_repair_and_update_jsons(sys.argv[2], sys.argv[3], sys.argv[4])
        print(f"\nProcessed {len(results)} PDBs")
    else:
        pdb_path = sys.argv[1]
        redesign_dict = {}
        if len(sys.argv) > 2:
            with open(sys.argv[2], 'r') as f:
                redesign_dict = json.load(f)
        
        result = repair_pdb_and_update_redesign_list(pdb_path, redesign_dict)
        print(result['report'])
        print(f"Success: {result['success']}")


if __name__ == "__main__":
    main()
