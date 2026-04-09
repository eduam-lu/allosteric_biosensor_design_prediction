"""

Given a backbone structure path, and a sequence or a set of sequences in a csv, 
it threads the sequence(s) onto the backbone and outputs the resulting structure(s) 
in the indicated output path

"""

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax

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
    
    return pose

# ==========================================
# Example Usage:
# ==========================================
if __name__ == "__main__":
    # Initialize PyRosetta
    pyrosetta.init("-mute all") # Muting output for cleaner console
    
    # Load your template backbone
    # For this example, replace 'template_backbone.pdb' with your actual file
    # template_pose = pyrosetta.pose_from_pdb("template_backbone.pdb")
    
    # Mock example (creating a pose from a sequence just to test the script)
    template_pose = pyrosetta.pose_from_sequence("ACDEFGHIKL") 
    
    # Define your new target sequence (must be the same length)
    my_new_sequence = "VWXYZVWXYZ" # Note: Z/X are not standard AAs, use valid 1-letter codes like "VWQYTVWQYT"
    my_new_sequence = "VWQYTVWQYT" 
    
    # Run the threading algorithm
    threaded_pose = thread_sequence(template_pose, my_new_sequence)