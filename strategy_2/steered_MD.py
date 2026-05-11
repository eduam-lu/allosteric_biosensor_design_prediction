import openmm as mm
import openmm.app as app
import openmm.unit as unit
import argparse
import os
import json

# ==========================================
# Parse Command-Line Arguments
# ==========================================
parser = argparse.ArgumentParser(description="Run Steered Molecular Dynamics simulation")
parser.add_argument('--input_pdb', required=True, help='Path to input equilibrated PDB file')
parser.add_argument('--output_dir', required=True, help='Directory to save output frames')
parser.add_argument('--domain_A_json', required=True, help='JSON file containing domain A atom indices')
parser.add_argument('--domain_B_json', required=True, help='JSON file containing domain B atom indices')
parser.add_argument('--total_time_ns', type=float, default=5.0, help='Total simulation time in nanoseconds (default: 5.0)')
parser.add_argument('--initial_distance', type=float, default=2.5, help='Initial COM distance in nm (default: 2.5)')
parser.add_argument('--target_distance', type=float, default=4.0, help='Target COM distance in nm (default: 4.0)')
parser.add_argument('--spring_constant', type=float, default=10000.0, help='Spring constant in kJ/mol/nm^2 (default: 10000.0)')
parser.add_argument('--num_frames', type=int, default=10, help='Number of frames to extract (default: 10)')
args = parser.parse_args()

# Create output directory
os.makedirs(args.output_dir, exist_ok=True)

# Load domain indices from JSON files
with open(args.domain_A_json, 'r') as f:
    domain_A_indices = json.load(f)
with open(args.domain_B_json, 'r') as f:
    domain_B_indices = json.load(f)

print(f"Domain A: {len(domain_A_indices)} atoms")
print(f"Domain B: {len(domain_B_indices)} atoms")

# ==========================================
# 1. System Setup & Parameters
# ==========================================
# Load pre-equilibrated structure (needs to be solvated/minimized first)
if not os.path.exists(args.input_pdb):
    raise FileNotFoundError(f"Input PDB not found: {args.input_pdb}")

pdb = app.PDBFile(args.input_pdb)

# Load standard Amber forcefield (adjust if using CHARMM)
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

system = forcefield.createSystem(pdb.topology, 
                                 nonbondedMethod=app.PME,
                                 nonbondedCutoff=1.0*unit.nanometer, 
                                 constraints=app.HBonds)

# ==========================================
# 2. Validate Domain Indices
# ==========================================
# Validate that indices are within bounds
total_atoms = pdb.topology.getNumAtoms()
if max(domain_A_indices) >= total_atoms or max(domain_B_indices) >= total_atoms:
    raise ValueError(f"Domain indices exceed total atoms ({total_atoms}) in PDB")
if min(domain_A_indices) < 0 or min(domain_B_indices) < 0:
    raise ValueError("Domain indices must be non-negative")

# ==========================================
# 3. Apply Steered MD Force (COM Distance)
# ==========================================
# CustomCentroidBondForce calculates the center of mass of groups and applies an equation
# k = spring constant, d0 = target distance (which we will change over time)
smd_equation = "0.5 * k * (distance(g1,g2) - d0)^2"
smd_force = mm.CustomCentroidBondForce(2, smd_equation)

# Define our variables as Global Parameters so we can update them mid-simulation
smd_force.addGlobalParameter("k", args.spring_constant) # kJ/mol/nm^2 (Stiff spring for pulling)
smd_force.addGlobalParameter("d0", args.initial_distance)    # nm (Starting distance)

# Add the groups to the force and link them
g1 = smd_force.addGroup(domain_A_indices)
g2 = smd_force.addGroup(domain_B_indices)
smd_force.addBond([g1, g2])
system.addForce(smd_force)

# ==========================================
# 4. Simulation Engine Setup
# ==========================================
# 300K, 1 ps^-1 friction, 2 fs timestep
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 
                                         1/unit.picosecond, 
                                         0.002*unit.picoseconds)

simulation = app.Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)

# Report progress to the console
simulation.reporters.append(app.StateDataReporter('stdout', 25000, step=True, 
                                                  potentialEnergy=True, temperature=True))

# ==========================================
# 5. The Steered Pulling Loop
# ==========================================
timestep_fs = 2.0
total_steps = int((args.total_time_ns * 1e6) / timestep_fs)

# We will update the distance every 1000 steps for a smooth pull
update_freq = 1000 
loops = total_steps // update_freq

# Calculate save frequency based on desired number of frames
save_freq_loops = max(1, loops // args.num_frames)

print(f"Starting {args.total_time_ns} ns SMD pull...")
print(f"Pulling from {args.initial_distance} nm to {args.target_distance} nm")
print(f"Spring constant: {args.spring_constant} kJ/mol/nm^2")
print(f"Total steps: {total_steps}")
print(f"Saving {args.num_frames} frames to {args.output_dir}")

for i in range(loops):
    # Calculate the new distance for this step
    fraction_complete = (i + 1) / loops
    current_d0 = args.initial_distance + (fraction_complete * (args.target_distance - args.initial_distance))
    
    # Update the parameter in the OpenMM Context dynamically
    simulation.context.setParameter("d0", current_d0)
    
    # Advance the simulation
    simulation.step(update_freq)
    
    # Extract the frame if we hit the target interval
    if (i + 1) % save_freq_loops == 0:
        frame_number = (i + 1) // save_freq_loops
        state = simulation.context.getState(getPositions=True)
        filename = os.path.join(args.output_dir, f'smd_frame_{frame_number}.pdb')
        with open(filename, 'w') as f:
            app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
        print(f"Saved frame {frame_number} at distance {current_d0:.2f} nm")

# print(f"\nSMD Pull Complete! Frames saved to {args.output_dir}")
# print("Ready for structure refinement and rescoring.")