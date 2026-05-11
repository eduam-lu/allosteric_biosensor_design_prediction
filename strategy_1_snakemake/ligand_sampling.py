"""
python ligand_sampling.py --ligand --structure --output --n_confs --n_pos --help

DESCRIPTION

This script, given an input ligand and a structure with a different ligand bound to it, will sample different conformations and positions for that new ligand to replace the other.
This is a previous step to binding site redesign with RF diffussion AA/ Ligand MPNN

FLAGS

--ligand: SMILES format (prolly a file is better)

--structure: path to pdb or cif with ligand bound to it

--output: where the output folders should be generated

--help: extended help

WORKFLOW

A. Sample conformers with RD kit and save them as pdb files

B. Identify ligand positions

1. Manually through moving the ligand in pymol and saving the coordinates as a pdb
2. Not manually, through alignment with the previous ligand and stochastic tilting, rotating and inverting 

C. Prepare final pdb files

1. Load a position
2. Load conformer
3. Align conformer with the position (maybe find a way to specify atoms specially important to align)
4. Iterate

D. Final visual adjustments if needed

Eduardo Amo González
2025-2026
"""

### IMPORT MODULES ##################################################################################

import argparse
import pandas as pd
from pathlib import Path
import sys
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
import numpy as np
import pymol

### PARAMS ##########################################################################################
num_conformers=5
conformer_rmsd_cutoff=0.75

num_positions=5
ligand_name="Tc"

rotation = 15
translation = 0.3
### FUNCTIONS #######################################################################################

def sample_conformers(molecule, n_conformers, rmsd_cutoff, output, ligand_name="lig"):
    output = Path(output)
    out_dir = output / "conformers"
    stats_file = out_dir / "conformer_stats.csv"
    
    # Check if the conformers actually FINISHED generating, not just if the folder exists
    if stats_file.exists():
        print("Conformers already found. Skipping conformer calculation.")
        df = pd.read_csv(stats_file)
        # Assuming sorted by energy or index, picking best
        best_idx = df.loc[df['energy'].idxmin()]['conf_id']
        return out_dir / f"{ligand_name}_conf_{int(best_idx)}.pdb"

    # Prepare output folder
    out_dir.mkdir(parents=True, exist_ok=True)

    # Load molecule
    mol = Chem.MolFromSmiles(molecule)
    if mol is None:
        raise ValueError(f"Could not parse SMILES: {molecule}")
    mol = Chem.AddHs(mol)

    # Generate conformers
    print("Generating conformers...")
    params = AllChem.ETKDG()
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)

    energies = []
    print("Minimizing conformers...")

    # Minimize each conformer & store energies
    for cid in conf_ids:
        prop = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94")
        ff = AllChem.MMFFGetMoleculeForceField(mol, prop, confId=cid)
        ff.Minimize()
        energy = ff.CalcEnergy()
        energies.append(energy)

    energies = np.array(energies)

    # Create heavy-atom-only copy for PDB output (RFD3 expects no hydrogens)
    mol_noH = Chem.RemoveAllHs(mol)

    # Set ligand residue name in PDB output
    for atom in mol_noH.GetAtoms():
        mi = Chem.AtomPDBResidueInfo()
        mi.SetResidueName(ligand_name[:3])
        mi.SetResidueNumber(1)
        mi.SetChainId("X")
        mi.SetIsHeteroAtom(True)
        symbol = atom.GetSymbol()
        idx = atom.GetIdx() + 1
        atom_name = f" {symbol}{idx}" if len(symbol) == 1 else f"{symbol}{idx}"
        mi.SetName(f"{atom_name:<4s}")
        atom.SetMonomerInfo(mi)

    # Save PDBs (heavy atoms only)
    for cid in conf_ids:
        filename = out_dir / f"{ligand_name}_conf_{cid}.pdb"
        Chem.MolToPDBFile(mol_noH, str(filename), confId=cid)

    print(f"{n_conformers} conformers generated and saved successfully.")

    # STATISTICS
    # Pairwise RMSD matrix
    rmsd_matrix = np.zeros((n_conformers, n_conformers))
    for a, i in enumerate(conf_ids):
        for b, j in enumerate(conf_ids):
            if b <= a: 
                continue
            rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=i, refId=j)
            rmsd_matrix[a, b] = rmsd
            rmsd_matrix[b, a] = rmsd

    # RMSD to lowest-energy conformer
    best_idx = energies.argmin()
    best_conf_id = conf_ids[int(best_idx)]
    lowest_energy_path = out_dir / f"{ligand_name}_conf_{best_conf_id}.pdb"

    rmsd_to_best = np.array([
        rdMolAlign.GetBestRMS(mol, mol, prbId=conf_ids[i], refId=best_conf_id)
        for i in range(n_conformers)
    ])

    # Count unique conformers using RMSD threshold
    unique = [0]  # always include lowest-energy conformer
    for i in range(1, n_conformers):
        if all(rmsd_matrix[i, j] > rmsd_cutoff for j in unique):
            unique.append(i)

    # Statistics output
    stats_df = pd.DataFrame({
        "conf_id": conf_ids,
        "energy": energies,
        "rmsd_to_lowest": rmsd_to_best,
        "is_unique": [cid in unique for cid in conf_ids]
    })

    stats_df.to_csv(out_dir / "conformer_stats.csv", index=False)
    np.savetxt(out_dir / "rmsd_matrix.txt", rmsd_matrix, fmt="%.3f")
    print(f"Unique conformers (RMSD > {rmsd_cutoff} Å): {len(unique)}")

    return lowest_energy_path



def identify_ligand_positions(molecule_path, structure, num_positions, output, rotation = rotation, translation=translation):
    output = Path(output)
    out_dir = output / "positions"
    
    if out_dir.exists() and any(out_dir.glob("position_*.pdb")):
        print("Positions already found. Skipping position calculation.")
        return

    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Limpieza inicial
    pymol.cmd.delete('all')

    # A. Cargar estructuras
    pymol.cmd.load(str(structure), 'receptor')
    pymol.cmd.load(str(molecule_path), 'new_ligand')

    # Identificar ligando viejo (dentro del receptor)
    pymol.cmd.select('old_ligand', 'receptor and organic')
    
    # Fallback si no hay 'organic'
    if pymol.cmd.count_atoms('old_ligand') == 0:
        pymol.cmd.select('old_ligand', 'receptor and hetatm and not resn HOH and not solvent')

    if pymol.cmd.count_atoms('old_ligand') == 0:
        print("Error: No ligand found in the reference structure.")
        return

    # B. ALINEAMIENTO ROBUSTO (CENTRO DE MASAS)
    # En lugar de 'fit' (que falla si los átomos son distintos), movemos el nuevo al centro del viejo.
    
    # 1. Obtener coordenadas del centro de masa (COM)
    com_old = pymol.cmd.centerofmass('old_ligand')
    com_new = pymol.cmd.centerofmass('new_ligand')
    
    # 2. Calcular vector de traslación (Destino - Origen)
    # Queremos mover NEW hacia OLD.
    translation_vector = [com_old[i] - com_new[i] for i in range(3)]
    
    # 3. Mover SOLO el nuevo ligando
    pymol.cmd.translate(translation_vector, 'new_ligand')
    
    print(f"Ligando inicial centrado en el sitio de unión (Shift: {translation_vector})")

    # C. Borrar el ligando viejo del objeto receptor para dejar el hueco libre
    pymol.cmd.remove('old_ligand')

    # D. Sampling suave
    print(f"Generating {num_positions} positions (Soft sampling)...")
    
    for i in range(num_positions):
        temp_lig_name = f'temp_lig_{i}'
        
        # Crear copia temporal del ligando ya centrado
        pymol.cmd.create(temp_lig_name, 'new_ligand')

        # --- PARÁMETROS SUAVES ---
        # Rotación: Máximo +/- 15 grados (perturbación ligera) en lugar de 360
        angle = np.random.uniform(-rotation, rotation)
        
        # Eje aleatorio
        axis = np.random.normal(size=3)
        if np.linalg.norm(axis) > 0: axis /= np.linalg.norm(axis)
        
        pymol.cmd.rotate(list(axis), angle, temp_lig_name)

        # Traslación: Máximo +/- 0.3 Angstroms (vibración local)
        trans = np.random.uniform(-translation, translation, size=3)
        pymol.cmd.translate(list(trans), temp_lig_name)
        
        # Guardar: Receptor (hueco) + Ligando Nuevo (movido)
        output_pdb = out_dir / f'position_{i}.pdb'
        pymol.cmd.save(str(output_pdb), f"receptor or {temp_lig_name}")
        
        # Limpieza iteración
        pymol.cmd.delete(temp_lig_name)

    # Limpieza final
    pymol.cmd.delete('all')
    print(f"{num_positions} positions were saved in {out_dir}")
    return

def insert_conformers(structure, output):
    output = Path(output)
    pos_dir = output / "positions"
    conf_dir = output / "conformers"
    final_dir = output / "final_pdbs"

    if not pos_dir.exists() or not conf_dir.exists():
        print("Missing positions or conformers directory.")
        return

    final_dir.mkdir(parents=True, exist_ok=True)
    
    positions = list(pos_dir.glob("position_*.pdb"))
    conformers = list(conf_dir.glob("*.pdb"))

    print(f"Procesando {len(positions)} posiciones x {len(conformers)} confórmeros...")

    for pos in positions:
        for conf in conformers:
            try:
                # 1. Limpiar sesión
                pymol.cmd.delete('all')

                # 2. Cargar estructuras
                # 'reference_complex': Contiene PROTEÍNA + LIGANDO (Posición)
                pymol.cmd.load(str(pos), 'reference_complex')
                
                # 'mobile_ligand': Contiene solo el LIGANDO (Conformero nuevo)
                pymol.cmd.load(str(conf), 'mobile_ligand')

                # 3. ALINEAMIENTO (FIT)
                # Ajustamos el 'mobile_ligand' sobre la parte 'organic' (ligando) del complejo de referencia.
                # 'organic' selecciona las moléculas pequeñas excluyendo la proteína.
                # Si tu ligando no es detectado como organic, usa 'hetatm'.
                try:
                    pymol.cmd.fit('mobile_ligand', 'reference_complex and organic', matchmaker=1)
                except pymol.CmdException:
                    print(f"Fit estricto falló para {conf.name}, intentando sin hidrógenos...")
                    # Intento secundario ignorando hidrógenos si hay discrepancias
                    pymol.cmd.fit('mobile_ligand and not hydro', 'reference_complex and organic and not hydro')

                # 4. PREPARAR ESTRUCTURA FINAL
                # Ya tenemos el 'mobile_ligand' en la posición correcta.
                # Ahora debemos borrar el ligando VIEJO de la proteína para que no choquen.
                pymol.cmd.remove('reference_complex and organic')

                # 5. CREAR EL NUEVO COMPLEJO
                # Unimos la proteína limpia (del reference_complex) + el nuevo ligando alineado (mobile_ligand)
                pymol.cmd.create('final_structure', 'reference_complex or mobile_ligand')

                # 6. Guardar
                pos_id = pos.stem.split('_')[-1]
                conf_id = conf.stem.split('_')[-1]
                output_pdb = final_dir / f"structure_pos{pos_id}_conf{conf_id}.pdb"
                
                pymol.cmd.save(str(output_pdb), 'final_structure')
                
            except Exception as e:
                print(f"Error procesando {pos.name} - {conf.name}: {e}")
                continue

    print(f"Final PDBs generated in {final_dir}")
    return

### MAIN PIPELINE FUNCTION ##########################################################################

def run_ligand_sampling_pipeline(ligand_smiles, structure_path, output_path, 
                          num_conformers=num_conformers, conformer_rmsd_cutoff=conformer_rmsd_cutoff, 
                          num_positions=num_positions, ligand_name=ligand_name, translation=translation, rotation=rotation):
    """
    Function to run the full pipeline programmatically.
    """
    structure_path = Path(structure_path)
    output_path = Path(output_path)

    # Validations
    if not structure_path.exists():
        raise FileNotFoundError(f"The file '{structure_path}' does not exist.")
    if structure_path.is_dir():
        raise IsADirectoryError(f"The path '{structure_path}' is a directory. File required")
    
    if not output_path.exists():
        output_path.mkdir(parents=True, exist_ok=True)

    print(f"Starting pipeline for ligand: {ligand_name}")
    print(f"Output directory: {output_path}")

    ### A. Sample conformers with RDkit
    lowest_energy_conformer = sample_conformers(
        ligand_smiles, 
        num_conformers, 
        conformer_rmsd_cutoff, 
        output_path, 
        ligand_name
    )

    ### B. Identify ligand positions
    identify_ligand_positions(
        lowest_energy_conformer, 
        structure_path, 
        num_positions, 
        output_path,
        rotation,
        translation,
    )

    ### C. Prepare final pdb files (insert conformer into each position)
    insert_conformers(
        structure_path,  
        output_path
    )
    
    print("Ligands sampled successfully.")
    return lowest_energy_conformer

### CLI EXECUTION ###################################################################################

if __name__ == "__main__":
    # This block only runs if the script is executed directly from terminal
    
    parser = argparse.ArgumentParser(
        description="This script runs the ligand sampling pipeline"
    )
    parser.add_argument('--ligand', help="SMILES code for the ligand", type=str)
    parser.add_argument('--structure', help="Path to the structure to redesign", type=str)
    parser.add_argument('--output', help="Folder where the outputs will be stored", type=str)
    # Optional arguments exposed to CLI
    parser.add_argument('--n_confs', help="Number of conformers", type=int, default=num_conformers)
    parser.add_argument('--n_pos', help="Number of positions", type=int, default=num_positions)

    parser.add_argument('--detailed-help', action='store_true', help="Show detailed help message and exit")

    args = parser.parse_args()

    # If --ligand, structure or output was not provided, show error and exit
    if args.detailed_help:
        parser.print_help()
        sys.exit(0)

    if not args.ligand or not args.structure or not args.output:
        parser.error("Arguments --ligand, --structure and --output are required.")

    # Run the pipeline using the function
    run_ligand_sampling_pipeline(
        ligand_smiles=args.ligand,
        structure_path=args.structure,
        output_path=args.output,
        num_conformers=args.n_confs,
        num_positions=args.n_pos,
        ligand_name="Tc", # Keep default or add CLI arg for this
        translation=translation,
        rotation=rotation,
        conformer_rmsd_cutoff=conformer_rmsd_cutoff
    )