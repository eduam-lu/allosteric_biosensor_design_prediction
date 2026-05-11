"""
This script takes a sequence, a backbone structure and a reference dimer and:
1. Threads the sequence onto the backbone structure
2. Aligns the threaded structure to the reference dimer
3. Extracts the interface residues from the aligned structure
4. Runs Rosetta's interface analyzer to assess the quality of the interface
5. Redesigns the interface residues using ProteinMPNN with a bias towards residues found in interfaces
6. If we are over the minimum iterations, and the interface quality has stopped improving or we have passed the number of iterations, stop. 
    Otherwise go back to step 1 with the redesigned sequence and the same backbone structure.
"""
### IMPORT MODULES #####################################################################
import os
import sys
import json
import argparse
import tempfile
import shutil
import subprocess
import numpy as np

import pyrosetta
from pyrosetta.rosetta.protocols.simple_moves import MutateResidue
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.analysis import InterfaceAnalyzerMover

### PARAMS and INPUT #####################################################################
parser = argparse.ArgumentParser(
    description="Iterative dimerization interface generator using LigandMPNN + Rosetta"
)
parser.add_argument("--sequence", type=str, required=True,
                    help="Input monomer amino acid sequence (1-letter code)")
parser.add_argument("--backbone", type=str, required=True,
                    help="Path to the monomer backbone PDB file")
parser.add_argument("--reference_dimer", type=str, required=True,
                    help="Path to the reference dimer PDB file (chains A and B)")
parser.add_argument("--output_folder", type=str, required=True,
                    help="Folder to write output structures and logs")
parser.add_argument("--max_iter", type=int, default=20,
                    help="Maximum number of optimisation iterations (default: 20)")
parser.add_argument("--min_iter", type=int, default=3,
                    help="Minimum iterations before early stopping is allowed (default: 3)")
parser.add_argument("--interface_cutoff", type=float, default=8.0,
                    help="Distance cutoff (Å) for defining interface residues (default: 8.0)")
parser.add_argument("--mpnn_temperature", type=float, default=0.1,
                    help="LigandMPNN sampling temperature (default: 0.1)")
parser.add_argument("--n_designs", type=int, default=10,
                    help="Number of LigandMPNN designs per iteration (default: 10)")
parser.add_argument("--nrelax", type=int, default=5,
                    help="Number of FastRelax trajectories (default: 5)")
parser.add_argument("--constraint_weight", type=float, default=1.0,
                    help="Coordinate constraint weight for FastRelax (default: 1.0)")
parser.add_argument("--reference_chain_a", type=str, default="A",
                    help="Chain ID for first monomer in reference dimer (default: A)")
parser.add_argument("--reference_chain_b", type=str, default="B",
                    help="Chain ID for second monomer in reference dimer (default: B)")
parser.add_argument("--ligandmpnn_path", type=str,
                    default=os.path.expanduser("~/LigandMPNN/run.py"),
                    help="Path to LigandMPNN run.py")
parser.add_argument("--ligandmpnn_checkpoint", type=str,
                    default=os.path.expanduser("~/LigandMPNN/model_params/ligandmpnn_v_32_010_25.pt"),
                    help="Path to LigandMPNN model checkpoint")
parser.add_argument("--ligandmpnn_conda_env", type=str,
                    default=os.path.expanduser("~/miniforge3/envs/ligandmpnn_env"),
                    help="Path to conda env for LigandMPNN")

args = parser.parse_args()

LIGANDMPNN_RUN   = args.ligandmpnn_path
LIGANDMPNN_CKPT  = args.ligandmpnn_checkpoint
LIGANDMPNN_ENV   = args.ligandmpnn_conda_env
### FUNCTIONS ######################################################################

def thread_sequence_onto_backbone(backbone_pose, target_sequence):
    """
    Thread a target sequence onto a backbone pose using MutateResidue + repack.
    Returns the threaded pose (modifies in place).
    """
    if len(target_sequence) != backbone_pose.total_residue():
        raise ValueError(
            f"Sequence length ({len(target_sequence)}) != pose length ({backbone_pose.total_residue()})"
        )
    for i, aa in enumerate(target_sequence):
        seqpos = i + 1
        if backbone_pose.sequence()[i] != aa:
            MutateResidue(seqpos, aa).apply(backbone_pose)

    scorefxn = pyrosetta.get_fa_scorefxn()
    task = pyrosetta.standard_packer_task(backbone_pose)
    task.restrict_to_repacking()
    PackRotamersMover(scorefxn, task).apply(backbone_pose)
    return backbone_pose


def align_to_reference(monomer_pose, reference_dimer_pose, ref_chain="A"):
    """
    Superimpose the monomer onto one chain of the reference dimer (CA atoms).
    Returns the aligned monomer pose.
    """
    from pyrosetta.rosetta.core.id import AtomID
    from pyrosetta.rosetta.protocols.toolbox.superimpose import superimpose_pose_on_subset_CA

    # Build a boolean vector: True for residues that belong to the target chain in the reference
    ref_info = reference_dimer_pose.pdb_info()
    ref_nres = reference_dimer_pose.total_residue()
    mono_nres = monomer_pose.total_residue()

    # Collect CA positions from the reference chain
    ref_ca_residues = []
    for r in range(1, ref_nres + 1):
        if ref_info.chain(r) == ref_chain:
            ref_ca_residues.append(r)

    if len(ref_ca_residues) != mono_nres:
        print(f"[WARN] Ref chain has {len(ref_ca_residues)} residues, monomer has {mono_nres}. "
              "Using min overlap for alignment.")

    n_align = min(len(ref_ca_residues), mono_nres)

    # Use pyrosetta's superimpose on subset
    from pyrosetta.rosetta.utility import vector1_bool
    atom_map = pyrosetta.rosetta.std.map_core_id_AtomID_core_id_AtomID()
    for i in range(n_align):
        mono_aid = AtomID(monomer_pose.residue(i + 1).atom_index("CA"), i + 1)
        ref_aid  = AtomID(reference_dimer_pose.residue(ref_ca_residues[i]).atom_index("CA"),
                          ref_ca_residues[i])
        atom_map[mono_aid] = ref_aid

    from pyrosetta.rosetta.core.scoring import superimpose_pose
    rmsd = superimpose_pose(monomer_pose, reference_dimer_pose, atom_map)
    print(f"  Alignment RMSD to reference chain {ref_chain}: {rmsd:.3f} Å")
    return monomer_pose


def build_dimer_pose(monomer_pose, reference_dimer_pose, ref_chain_b="B"):
    """
    Create a dimer by combining the aligned monomer (as chain A) with chain B 
    extracted from the reference dimer. Returns the dimer pose.
    """
    ref_info = reference_dimer_pose.pdb_info()
    ref_nres = reference_dimer_pose.total_residue()

    # Extract chain B residues from reference
    chain_b_start = None
    chain_b_end = None
    for r in range(1, ref_nres + 1):
        if ref_info.chain(r) == ref_chain_b:
            if chain_b_start is None:
                chain_b_start = r
            chain_b_end = r

    if chain_b_start is None:
        raise ValueError(f"Chain {ref_chain_b} not found in reference dimer")

    # Clone monomer as chain A, append chain B from reference
    dimer_pose = monomer_pose.clone()
    from pyrosetta.rosetta.core.pose import append_subpose_to_pose
    append_subpose_to_pose(dimer_pose, reference_dimer_pose, chain_b_start, chain_b_end, True)

    return dimer_pose


def get_interface_residues(dimer_pose, cutoff=8.0):
    """
    Identify interface residues (on chain A side) within `cutoff` Å of chain B.
    Returns a list of residue numbers (1-indexed, chain A).
    """
    info = dimer_pose.pdb_info()
    nres = dimer_pose.total_residue()

    chain_a_res = [r for r in range(1, nres + 1) if info.chain(r) == "A"]
    chain_b_res = [r for r in range(1, nres + 1) if info.chain(r) == "B"]

    interface_res = []
    for ra in chain_a_res:
        ca_a = dimer_pose.residue(ra).xyz("CA")
        for rb in chain_b_res:
            ca_b = dimer_pose.residue(rb).xyz("CA")
            dist = ca_a.distance(ca_b)
            if dist <= cutoff:
                interface_res.append(ra)
                break

    return interface_res


def interface_analyzer(dimer_pose, interface="A_B"):
    """
    Run Rosetta InterfaceAnalyzerMover on the dimer and return interface metrics.
    """
    scorefxn = pyrosetta.get_fa_scorefxn()
    iam = InterfaceAnalyzerMover(interface)
    iam.set_scorefunction(scorefxn)
    iam.set_compute_packstat(True)
    iam.set_pack_separated(True)
    iam.apply(dimer_pose)

    dG = dimer_pose.scores.get("dG_separated", float("inf"))
    dSASA = dimer_pose.scores.get("dSASA_int", 0.0)
    packstat = dimer_pose.scores.get("packstat", 0.0)

    print(f"  Interface dG_separated: {dG:.2f} | dSASA: {dSASA:.1f} | packstat: {packstat:.3f}")
    return {"dG_separated": dG, "dSASA_int": dSASA, "packstat": packstat}


def relax_dimer(dimer_pose, constraint_weight=1.0, nsamples=5):
    """
    Run FastRelax on the dimer with coordinate constraints.
    Returns the best-scoring relaxed pose.
    """
    scorefxn = pyrosetta.get_fa_scorefxn()
    scorefxn.set_weight(
        pyrosetta.rosetta.core.scoring.coordinate_constraint, constraint_weight
    )

    best_pose = None
    best_score = float("inf")

    for i in range(nsamples):
        trial_pose = dimer_pose.clone()
        relax = FastRelax(scorefxn, 15)
        relax.constrain_relax_to_start_coords(True)
        relax.apply(trial_pose)
        score = scorefxn(trial_pose)
        if score < best_score:
            best_score = score
            best_pose = trial_pose.clone()
        print(f"    Relax sample {i+1}/{nsamples}: score = {score:.2f}")

    print(f"  Best relax score: {best_score:.2f}")
    return best_pose


def run_ligandmpnn_biased(pdb_path, interface_residues, output_folder,
                          temperature=0.1, n_designs=10, chain_id="A"):
    """
    Run LigandMPNN with redesign restricted to interface residues.
    Returns a list of generated sequence strings.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Build redesigned_residues string: "A5 A6 A12 ..."
    res_str = " ".join(f"{chain_id}{r}" for r in interface_residues)

    # Write JSON inputs
    paths_json = {pdb_path: ""}
    residues_json = {pdb_path: res_str}

    paths_file = os.path.join(output_folder, "pdb_paths.json")
    res_file = os.path.join(output_folder, "redesigned_residues.json")
    with open(paths_file, "w") as f:
        json.dump(paths_json, f)
    with open(res_file, "w") as f:
        json.dump(residues_json, f)

    cmd = [
        "conda", "run", "-p", LIGANDMPNN_ENV, "--no-banner",
        "python", LIGANDMPNN_RUN,
        "--pdb_path_multi", paths_file,
        "--out_folder", output_folder,
        "--save_stats", "1",
        "--batch_size", str(n_designs),
        "--number_of_batches", "1",
        "--model_type", "ligand_mpnn",
        "--checkpoint_ligand_mpnn", LIGANDMPNN_CKPT,
        "--temperature", str(temperature),
        "--redesigned_residues_multi", res_file,
    ]

    print(f"  Running LigandMPNN on {len(interface_residues)} interface residues...")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  [ERROR] LigandMPNN failed:\n{result.stderr}")
        return []

    # Parse output FASTA files
    seqs_dir = os.path.join(output_folder, "seqs")
    sequences = []
    if os.path.isdir(seqs_dir):
        for fasta_file in sorted(os.listdir(seqs_dir)):
            if fasta_file.endswith(".fa"):
                with open(os.path.join(seqs_dir, fasta_file)) as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith(">"):
                            sequences.append(line)

    print(f"  LigandMPNN generated {len(sequences)} sequences")
    return sequences


def pick_best_sequence(sequences, backbone_pose, reference_dimer_pose,
                       ref_chain_a, ref_chain_b, interface_cutoff, constraint_weight):
    """
    Thread each candidate sequence, build dimer, quick-relax, score interface.
    Returns (best_sequence, best_dG, best_pose).
    """
    best_seq = None
    best_dG = float("inf")
    best_pose = None

    for idx, seq in enumerate(sequences):
        print(f"  Evaluating candidate {idx+1}/{len(sequences)}...")
        trial_pose = backbone_pose.clone()
        try:
            thread_sequence_onto_backbone(trial_pose, seq)
        except ValueError as e:
            print(f"    Skipping: {e}")
            continue

        align_to_reference(trial_pose, reference_dimer_pose, ref_chain=ref_chain_a)
        dimer = build_dimer_pose(trial_pose, reference_dimer_pose, ref_chain_b=ref_chain_b)

        # Quick single-trajectory relax for screening
        scorefxn = pyrosetta.get_fa_scorefxn()
        scorefxn.set_weight(
            pyrosetta.rosetta.core.scoring.coordinate_constraint, constraint_weight
        )
        relax = FastRelax(scorefxn, 5)
        relax.constrain_relax_to_start_coords(True)
        relax.apply(dimer)

        metrics = interface_analyzer(dimer)
        if metrics["dG_separated"] < best_dG:
            best_dG = metrics["dG_separated"]
            best_seq = seq
            best_pose = dimer.clone()

    return best_seq, best_dG, best_pose


### MAIN EXECUTION ######################################################################
def main():
    pyrosetta.init("-ignore_unrecognized_res -load_PDB_components false -mute all")

    os.makedirs(args.output_folder, exist_ok=True)

    # Load structures
    print("Loading structures...")
    backbone_pose = pyrosetta.pose_from_file(args.backbone)
    reference_dimer_pose = pyrosetta.pose_from_file(args.reference_dimer)

    # Step 1: Thread initial sequence onto backbone
    print("\n=== Iteration 0: Initial threading ===")
    current_sequence = args.sequence
    current_pose = backbone_pose.clone()
    thread_sequence_onto_backbone(current_pose, current_sequence)

    # Step 2: Align to reference dimer
    align_to_reference(current_pose, reference_dimer_pose, ref_chain=args.reference_chain_a)

    # Step 3: Build initial dimer
    dimer_pose = build_dimer_pose(current_pose, reference_dimer_pose,
                                  ref_chain_b=args.reference_chain_b)

    # Step 4: Relax dimer
    print("  Relaxing initial dimer...")
    dimer_pose = relax_dimer(dimer_pose, constraint_weight=args.constraint_weight,
                             nsamples=args.nrelax)

    # Step 5: Score initial interface
    best_metrics = interface_analyzer(dimer_pose)
    best_dG = best_metrics["dG_separated"]
    best_sequence = current_sequence
    best_dimer_pose = dimer_pose.clone()

    # Save initial dimer
    init_path = os.path.join(args.output_folder, "iter_0_dimer.pdb")
    dimer_pose.dump_pdb(init_path)
    print(f"  Initial dimer saved to {init_path}")

    # Iterative optimisation loop
    no_improvement_count = 0
    results_log = [{"iteration": 0, "dG_separated": best_dG, "sequence": best_sequence}]

    for iteration in range(1, args.max_iter + 1):
        print(f"\n=== Iteration {iteration}/{args.max_iter} ===")

        # Identify interface residues from current best dimer
        interface_res = get_interface_residues(best_dimer_pose, cutoff=args.interface_cutoff)
        if not interface_res:
            print("  No interface residues found. Stopping.")
            break
        print(f"  Interface residues ({len(interface_res)}): {interface_res[:10]}{'...' if len(interface_res) > 10 else ''}")

        # Run biased LigandMPNN on interface residues
        mpnn_out = os.path.join(args.output_folder, f"mpnn_iter_{iteration}")
        dimer_pdb = os.path.join(args.output_folder, f"iter_{iteration-1}_dimer.pdb")
        if not os.path.exists(dimer_pdb):
            best_dimer_pose.dump_pdb(dimer_pdb)

        candidate_sequences = run_ligandmpnn_biased(
            dimer_pdb, interface_res, mpnn_out,
            temperature=args.mpnn_temperature, n_designs=args.n_designs
        )

        if not candidate_sequences:
            print("  No sequences generated. Continuing with current best.")
            no_improvement_count += 1
            if no_improvement_count >= 3 and iteration >= args.min_iter:
                print("  Stopping: 3 consecutive iterations without new candidates.")
                break
            continue

        # Evaluate candidates
        new_seq, new_dG, new_dimer = pick_best_sequence(
            candidate_sequences, backbone_pose, reference_dimer_pose,
            args.reference_chain_a, args.reference_chain_b,
            args.interface_cutoff, args.constraint_weight
        )

        if new_seq is None:
            print("  All candidates failed evaluation.")
            no_improvement_count += 1
        elif new_dG < best_dG:
            print(f"  IMPROVEMENT: dG {best_dG:.2f} -> {new_dG:.2f}")
            best_dG = new_dG
            best_sequence = new_seq
            no_improvement_count = 0

            # Full relax on best candidate
            best_dimer_pose = relax_dimer(new_dimer, constraint_weight=args.constraint_weight,
                                          nsamples=args.nrelax)
            # Re-score after full relax
            best_metrics = interface_analyzer(best_dimer_pose)
            best_dG = best_metrics["dG_separated"]

            out_path = os.path.join(args.output_folder, f"iter_{iteration}_dimer.pdb")
            best_dimer_pose.dump_pdb(out_path)
            print(f"  Saved improved dimer to {out_path}")
        else:
            print(f"  No improvement: best dG = {best_dG:.2f}, candidate dG = {new_dG:.2f}")
            no_improvement_count += 1

        results_log.append({
            "iteration": iteration,
            "dG_separated": best_dG,
            "sequence": best_sequence,
        })

        # Early stopping check
        if no_improvement_count >= 3 and iteration >= args.min_iter:
            print(f"\n  Early stopping after {no_improvement_count} iterations without improvement.")
            break

    # Save final outputs
    final_pdb = os.path.join(args.output_folder, "best_dimer.pdb")
    best_dimer_pose.dump_pdb(final_pdb)

    final_log = os.path.join(args.output_folder, "optimisation_log.json")
    with open(final_log, "w") as f:
        json.dump(results_log, f, indent=2)

    final_seq_file = os.path.join(args.output_folder, "best_sequence.txt")
    with open(final_seq_file, "w") as f:
        f.write(best_sequence + "\n")

    print(f"\n{'='*60}")
    print(f"DONE — Best dG_separated: {best_dG:.2f}")
    print(f"  Best dimer:    {final_pdb}")
    print(f"  Best sequence: {final_seq_file}")
    print(f"  Log:           {final_log}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()