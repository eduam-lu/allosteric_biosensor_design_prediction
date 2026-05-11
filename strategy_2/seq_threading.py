"""

This script accepts a csv file with two columns: "backbone_path" and "target_sequence". 
It reads each backbone structure, threads the corresponding target sequence onto it, and outputs the threaded structure as a PDB file. 
The output filename is derived from the input backbone filename and the target sequence for clarity.

The script returns a csv file with the original "backbone_path", "target_sequence", and thenew column "threaded_pdb_path" that contains the path to the newly created threaded PDB file
and "energy" that contains the score of the threaded structure after packing and relaxation.
[USAGE]
seq_threading.py --input_csv input_sequences.csv --output_csv threaded_output.csv --output_dir threaded_structures/ --num_workers 8

"""

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
import argparse
import csv
import os
from pathlib import Path
from multiprocessing import Pool, cpu_count

def thread_sequence(pose, target_sequence, output_filename="threaded_output.pdb"):
    """
    Threads a target sequence onto a given backbone pose.
    """
    # 1. Verify lengths match
    if len(target_sequence) != pose.total_residue():
        raise ValueError(f"Sequence length ({len(target_sequence)}) does not match Pose length ({pose.total_residue()}).")
    
    print(f"Threading sequence of length {len(target_sequence)}...")

    # 2. Mutate each residue to the target sequence
    for i, aa in enumerate(target_sequence):
        seqpos = i + 1  # PyRosetta is 1-indexed
        
        # Only mutate if the amino acid is different to save computation
        if pose.sequence()[i] != aa:
            mutator = MutateResidue(seqpos, aa)
            mutator.apply(pose)
            
    print("Mutations complete. Initializing packing and relaxation...")

    # Set up the standard full-atom score function
    scorefxn = pyrosetta.get_fa_scorefxn()

    # 3. Repack Sidechains
    # Create a task factory and restrict it to repacking only (no further mutations)
    task_pack = pyrosetta.standard_packer_task(pose)
    task_pack.restrict_to_repacking()
    
    pack_mover = PackRotamersMover(scorefxn, task_pack)
    pack_mover.apply(pose)
    print("Repacking complete.")

    # 4. Relax the structure to resolve remaining backbone/sidechain clashes
    relax = FastRelax()
    relax.set_scorefxn(scorefxn)
    
    # Optional: You can constrain the backbone to its original coordinates 
    # if you want to strictly maintain the exact template backbone.
    relax.constrain_relax_to_start_coords(True)
    
    relax.apply(pose)
    print("Relaxation complete.")

    # 5. Output the final threaded structure
    pose.dump_pdb(output_filename)
    print(f"Threaded structure saved to {output_filename}")
    
    # 6. Calculate and return the final energy score
    energy = scorefxn(pose)
    print(f"Final energy score: {energy}")
    
    return pose, energy

# ==========================================
# Worker Function for Parallel Processing:
# ==========================================
def thread_single_sequence(backbone_path, target_sequence, output_dir):
    """
    Worker function that processes a single sequence in a separate process.
    PyRosetta is initialized fresh in each worker.
    """
    import pyrosetta
    from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
    from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pathlib import Path
    import os
    
    # Initialize PyRosetta in this worker process
    pyrosetta.init("-mute all")
    
    try:
        if not os.path.exists(backbone_path):
            return {
                'backbone_path': backbone_path,
                'target_sequence': target_sequence,
                'threaded_pdb_path': 'ERROR',
                'energy': 'ERROR',
                'error': f'Backbone not found: {backbone_path}'
            }
        
        pose = pyrosetta.pose_from_pdb(backbone_path)
        
        if len(target_sequence) != pose.total_residue():
            return {
                'backbone_path': backbone_path,
                'target_sequence': target_sequence,
                'threaded_pdb_path': 'ERROR',
                'energy': 'ERROR',
                'error': f'Sequence length mismatch'
            }
        
        # Perform mutations
        for i, aa in enumerate(target_sequence):
            seqpos = i + 1
            if pose.sequence()[i] != aa:
                mutator = MutateResidue(seqpos, aa)
                mutator.apply(pose)
        
        # Packing and relaxation
        scorefxn = pyrosetta.get_fa_scorefxn()
        task_pack = pyrosetta.standard_packer_task(pose)
        task_pack.restrict_to_repacking()
        pack_mover = PackRotamersMover(scorefxn, task_pack)
        pack_mover.apply(pose)
        
        relax = FastRelax()
        relax.set_scorefxn(scorefxn)
        relax.constrain_relax_to_start_coords(True)
        relax.apply(pose)
        
        # Save output
        backbone_name = Path(backbone_path).stem
        output_filename = os.path.join(output_dir, f"{backbone_name}_{target_sequence}.pdb")
        pose.dump_pdb(output_filename)
        
        energy = scorefxn(pose)
        
        return {
            'backbone_path': backbone_path,
            'target_sequence': target_sequence,
            'threaded_pdb_path': output_filename,
            'energy': energy
        }
    
    except Exception as e:
        return {
            'backbone_path': backbone_path,
            'target_sequence': target_sequence,
            'threaded_pdb_path': 'ERROR',
            'energy': 'ERROR',
            'error': str(e)
        }

# ==========================================
# Main Processing Function (Parallel):
# ==========================================
def process_sequences_from_csv(input_csv, output_csv, output_dir, num_workers=None):
    """
    Reads sequence threading jobs from input CSV, processes them in parallel, and writes results to output CSV.
    
    Args:
        input_csv (str): Path to input CSV with columns "backbone_path" and "target_sequence"
        output_csv (str): Path to output CSV to write results
        output_dir (str): Directory to save threaded PDB structures
        num_workers (int): Number of parallel workers (default: all available CPUs)
    """
    if num_workers is None:
        num_workers = cpu_count()
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Read input CSV
    input_data = []
    try:
        with open(input_csv, 'r') as f:
            reader = csv.DictReader(f)
            input_data = list(reader)
        print(f"Loaded {len(input_data)} sequences from {input_csv}")
    except FileNotFoundError:
        print(f"Error: Input CSV file '{input_csv}' not found.")
        return
    except Exception as e:
        print(f"Error reading input CSV: {e}")
        return
    
    # Prepare tasks (filter out empty rows)
    tasks = []
    for row in input_data:
        backbone_path = row.get('backbone_path', '').strip()
        target_sequence = row.get('target_sequence', '').strip()
        
        if backbone_path and target_sequence:
            tasks.append((backbone_path, target_sequence, output_dir))
    
    if not tasks:
        print("No valid tasks found in input CSV.")
        return
    
    # Process in parallel
    print(f"Starting parallel processing with {num_workers} workers...")
    with Pool(num_workers) as pool:
        results = pool.starmap(thread_single_sequence, tasks)
    
    # Write output CSV
    try:
        with open(output_csv, 'w', newline='') as f:
            fieldnames = ['backbone_path', 'target_sequence', 'threaded_pdb_path', 'energy']
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(results)
        print(f"\nResults written to {output_csv}")
        print(f"Completed {len(results)} sequences")
    except Exception as e:
        print(f"Error writing output CSV: {e}")


# ==========================================
# CLI Argument Parsing and Example Usage:
# ==========================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Thread target sequences onto backbone structures using PyRosetta (parallel processing)."
    )
    parser.add_argument(
        '--input_csv',
        required=True,
        help='Path to input CSV file with columns "backbone_path" and "target_sequence"'
    )
    parser.add_argument(
        '--output_csv',
        required=True,
        help='Path to output CSV file with results'
    )
    parser.add_argument(
        '--output_dir',
        required=True,
        help='Directory to save threaded PDB structures'
    )
    parser.add_argument(
        '--num_workers',
        type=int,
        default=None,
        help='Number of parallel workers (default: all available CPUs)'
    )
    
    args = parser.parse_args()
    
    # Process sequences in parallel
    process_sequences_from_csv(args.input_csv, args.output_csv, args.output_dir, num_workers=args.num_workers)