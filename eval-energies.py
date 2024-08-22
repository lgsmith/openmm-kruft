from openmm.app import Simulation, PDBFile
import openmm as mm
import mdtraj as md
from openmm import unit as u


import argparse as ap
import json
import numpy as np
from pathlib import Path

p = ap.ArgumentParser(formatter_class=ap.ArgumentDefaultsHelpFormatter)
p.add_argument('system', type=Path,
               help='system.xml to create an openmm simulation from.')
p.add_argument('dcd', type=str,
               help='Trajectory to re-evaluate energies along.')
p.add_argument('pdb_with_bonds', type=str,
               help='Path to a topology-specifying pdb with which the traj can be read.'
               ' Should also match the system.xml in atom order.')
p.add_argument('potential_energies', type=Path, 
               help='Name of file to write energies to. If .txt provided, will be a text file. Otherwise will save in numpy binary format.')
p.add_argument('--with-amber14-tip3p', action=ap.BooleanOptionalAction,
               help='If thrown, build system from pdb given as "topology"'
               ' using tip3p and amber14 from openmm.')
p.add_argument('--platform', type=str, default=None,
               help='If provided, use specified platform.')
p.add_argument('--platform-properties', type=Path, default=None,
               help='If provided, a json text file containing the desired energy properties.')

args = p.parse_args()
# deserialize system and open traj
sys = mm.XmlSerializer.deserialize(args.system.read_text())
top = PDBFile(args.pdb_with_bonds)
traj = md.load_dcd(args.dcd, top=args.pdb_with_bonds)
teg = mm.VerletIntegrator(0.002)
# Figure out if we're doing a custom platform here:
if args.platform:
    platform = mm.Platform.getPlatformByName(args.platform)
else:
    platform = None
if args.platform_properties:
    platform_properties = json.loads(args.platform_properties.read_text())
else:
    platform_properties = None

# set up simulation
sim = Simulation(top, sys, teg, platform=args.platform, platformProperties=platform_properties)
print('number of particles in Simulation:', sim.context.getSystem().getNumParticles())
print('shape of traj.xyz', traj.xyz.shape)
# set up empty numpy array; fill with energies from traj
potential_energies = np.zeros(len(traj))
for i, frame_crds in enumerate(traj.xyz):
    sim.context.setPositions(frame_crds)
    # if you also wanted to look at forces you'd request that here.
    state = sim.context.getState(getEnergy=True)
    ene = state.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
    print('frame:', 'energy (kcal/mol)', ene)
    potential_energies[i] = ene
# write energies to file
pe_out = args.potential_energies
if pe_out.suffix == '.txt':
    np.savetxt(pe_out, potential_energies)
else:
    if pe_out.suffix != '.npy':
        print('Warning, this script is going to save a numpy binary file to', 
              pe_out, 
              'but you have chosen an extension that is not ".npy". Proceeding.')
    np.save(pe_out, potential_energies)

