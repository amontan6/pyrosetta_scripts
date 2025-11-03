import sys
import argparse
from pyrosetta import *
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import MonteCarlo

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser(description="PyRosetta file loader example")

parser.add_argument("filename", help="Path to the input PDB file")
args = parser.parse_args()

print(f"Filename: {args.filename}")

mypose = pose_from_pdb(args.filename)

print(f"Loaded pose with {mypose.total_residue()} residues from: {args.filename}")

scorefxn = get_score_function(True)
score = scorefxn(mypose)
print(score)

temperature = 1.0
mc = MonteCarlo(mypose, scorefxn, temperature) ##come back and check this sytax

for i in range(1):
    n_res = mypose.total_residue()
    random_res = numeric.random.uniform()
    residue_index = int(random_res * n_res) + 1

    phi_perturb = numeric.random.gaussian()
    psi_perturb = numeric.random.gaussian()
    orig_phi = mypose.phi(residue_index)
    orig_psi = mypose.psi(residue_index)
    mypose.set_phi(residue_index, orig_phi + phi_perturb)
    mypose.set_psi(residue_index, orig_psi + psi_perturb)

    accepted = mc.boltzmann(mypose)
    print(accepted)

