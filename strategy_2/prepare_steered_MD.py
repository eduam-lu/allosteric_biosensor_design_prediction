import openmm as mm
import openmm.app as app
import openmm.unit as unit
import argparse
import os

# Parse command-line arguments
parser = argparse.ArgumentParser(description="Prepare protein structure for Steered MD simulation")
parser.add_argument('--input_pdb', required=True, help='Path to input PDB file')
parser.add_argument('--output_pdb', default='equilibrated_system.pdb', help='Path to output PDB file (default: equilibrated_system.pdb)')
args = parser.parse_args()

# Verify input file exists
if not os.path.exists(args.input_pdb):
    raise FileNotFoundError(f"Input PDB file not found: {args.input_pdb}")

# Load input PDB
pdb = app.PDBFile(args.input_pdb)

# Standard Amber14 forcefield and TIP3P water model
forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Modeller object 
modeller = app.Modeller(pdb.topology, pdb.positions)

#  Add Hydrogens at pH 7.0...
modeller.addHydrogens(forcefield, pH=7.0)

# Building the Water Box and Adding Ions 
# padding=1.0*nm means the box walls are 1 nm away from the protein on all sides.
# ionicStrength=0.15*molar mimics physiological salt (NaCl) and automatically neutralizes the system.

modeller.addSolvent(forcefield, padding=1.0*unit.nanometers, ionicStrength=0.15*unit.molar)

# Creating the System object

# Now we convert the topology into a physical math model
system = forcefield.createSystem(modeller.topology, 
                                 nonbondedMethod=app.PME, 
                                 nonbondedCutoff=1.0*unit.nanometer, 
                                 constraints=app.HBonds)

# Set up the integrator (Standard 300K simulation)
integrator = mm.LangevinMiddleIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)

# Use GPU if available, otherwise fallback to CPU
try:
    platform = mm.Platform.getPlatformByName('CUDA')
    properties = {'Precision': 'mixed'}
    simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
except Exception:
    print("CUDA not found, falling back to CPU...")
    simulation = app.Simulation(modeller.topology, system, integrator)

simulation.context.setPositions(modeller.positions)

# Minimizing Energy 

# This is needed in case the rosetta threading wasn't perfect

simulation.minimizeEnergy()

# Equilibrating the Water (100 ps)

# Let the water molecules bump around and settle into the protein crevices
# 50,000 steps at 2fs = 100 picoseconds

simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
simulation.step(50000) 

# Saving equilibrated_system.pdb

# Extract the relaxed state
state = simulation.context.getState(getPositions=True)

with open(args.output_pdb, 'w') as f:
    app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

#print(f"System prepared and saved to {args.output_pdb}")