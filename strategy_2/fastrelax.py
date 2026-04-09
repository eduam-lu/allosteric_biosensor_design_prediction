"""
Given a structure path, it returns a minimised structure in the indicated output path.
Runs N parallel FastRelax trajectories and keeps the best-scoring result.
"""
### IMPORT MODULES
import os
import csv
import argparse
import tempfile
import shutil
from multiprocessing import Pool

### ARGPARSE SETUP
parser = argparse.ArgumentParser(
    description="Relax a PDB/CIF structure or all structures in a folder. Outputs relaxed structures to a given folder."
)
parser.add_argument('--input_structure', help="Path to input .cif/.pdb file or folder with structure files", type=str, required=True)
parser.add_argument('--output_folder', help="Folder to output the relaxed structures", type=str, required=True)
parser.add_argument('--constraint_weight', help="Weight for coordinate constraints (default: 1.0)", type=float, default=1.0)
parser.add_argument('--nsamples', help="Number of parallel relax trajectories per structure (default: 5)", type=int, default=5)
parser.add_argument('--ncpus', help="Number of CPUs for parallelisation (default: all available)", type=int, default=None)

args = parser.parse_args()
input_path = args.input_structure
output_folder = args.output_folder
constraint_weight = args.constraint_weight
nsamples = args.nsamples
ncpus = args.ncpus

### WORKER INITIALIZER — each subprocess gets its own PyRosetta instance
def _init_worker():
    import pyrosetta
    pyrosetta.init("-ignore_unrecognized_res -load_PDB_components false -mute all")

### SINGLE RELAX TRAJECTORY (runs in a subprocess)
def _relax_one(args_tuple):
    pdb_path, sample_idx, tmp_dir, cst_weight = args_tuple
    import pyrosetta
    from pyrosetta.rosetta.protocols.relax import FastRelax
    from pyrosetta.rosetta.protocols.idealize import IdealizeMover

    file_name = os.path.basename(pdb_path)
    base_name = os.path.splitext(file_name)[0]
    ext = os.path.splitext(file_name)[1].lower()
    tmp_path = os.path.join(tmp_dir, f"{base_name}_sample{sample_idx}{ext}")

    try:
        pose = pyrosetta.pose_from_file(pdb_path)
        IdealizeMover().apply(pose)
        scorefxn = pyrosetta.get_fa_scorefxn()
        scorefxn.set_weight(pyrosetta.rosetta.core.scoring.coordinate_constraint, cst_weight)
        relax = FastRelax(scorefxn, 15)
        relax.constrain_relax_to_start_coords(True)
        relax.apply(pose)
        score = scorefxn(pose)
        if ext == '.cif':
            pose.dump_cif(tmp_path)
        else:
            pose.dump_pdb(tmp_path)
        return (score, tmp_path)
    except Exception as e:
        print(f"[!] Sample {sample_idx} failed on {file_name}: {e}")
        return (float('inf'), None)

### RELAX SINGLE FILE WITH N SAMPLES
def relax_single(pdb_path, output_dir, cst_weight=1.0, n=5, cpus=None):
    os.makedirs(output_dir, exist_ok=True)
    file_name = os.path.basename(pdb_path)
    base_name = os.path.splitext(file_name)[0]
    ext = os.path.splitext(file_name)[1].lower()
    output_path = os.path.join(output_dir, f"{base_name}_relaxed{ext}")

    tmp_dir = tempfile.mkdtemp(prefix=f"relax_{base_name}_")
    worker_args = [(pdb_path, i, tmp_dir, cst_weight) for i in range(n)]

    workers = min(cpus or os.cpu_count(), n)
    with Pool(processes=workers, initializer=_init_worker) as pool:
        results = pool.map(_relax_one, worker_args)

    # Select best scoring result
    results = [(s, p) for s, p in results if p is not None]
    if not results:
        print(f"[!] All {n} samples failed for {file_name}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        return None

    best_score, best_path = min(results, key=lambda x: x[0])
    shutil.copy2(best_path, output_path)
    shutil.rmtree(tmp_dir, ignore_errors=True)
    print(f"[OK] {file_name} -> best score: {best_score:.2f} (from {n} samples)")
    return best_score

### FUNCTION TO HANDLE DIRECTORY
def relax_folder(input_dir, output_dir, cst_weight=1.0, n=5, cpus=None):
    struct_files = sorted(f for f in os.listdir(input_dir) if f.lower().endswith(('.cif','.pdb')))
    if not struct_files:
        print("[!] No .cif or .pdb files found.")
        return

    scores = {}
    for struct in struct_files:
        full_path = os.path.join(input_dir, struct)
        best = relax_single(full_path, output_dir, cst_weight, n, cpus)
        if best is not None:
            scores[struct] = best

    # Write summary CSV
    csv_path = os.path.join(output_dir, "relax_scores.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["structure", "best_score"])
        for name, score in scores.items():
            writer.writerow([name, f"{score:.4f}"])
    print(f"[OK] Scores written to {csv_path}")
    return scores

### MAIN LOGIC
if os.path.isfile(input_path) and input_path.lower().endswith(('.cif', '.pdb')):
    best = relax_single(input_path, output_folder, constraint_weight, nsamples, ncpus)
    if best is not None:
        csv_path = os.path.join(output_folder, "relax_scores.csv")
        with open(csv_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["structure", "best_score"])
            writer.writerow([os.path.basename(input_path), f"{best:.4f}"])
        print(f"[OK] Score written to {csv_path}")
elif os.path.isdir(input_path):
    relax_folder(input_path, output_folder, constraint_weight, nsamples, ncpus)
else:
    print(f"[!] Input must be a .cif/.pdb file or a folder containing structure files. Got: {input_path}")