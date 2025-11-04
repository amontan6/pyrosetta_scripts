import sys
import argparse
from pyrosetta import *
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.kinematics import MoveMap
from pyrosetta.rosetta.core.optimization import MinimizerOptions, AtomTreeMinimizer
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers


init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser(description="PyRosetta file loader example")

parser.add_argument("filename", help="Path to the input PDB file")
args = parser.parse_args()

print(f"Filename: {args.filename}")

mypose = pose_from_pdb(args.filename)

print(f"Loaded pose with {mypose.total_residue()} residues from: {args.filename}")

scorefxn = get_score_function(True)
current_score = scorefxn(mypose)
print(current_score)
print("-------------")

temperature = 1.0
mc = MonteCarlo(mypose, scorefxn, temperature) ##come back and check this syntax

# Set up MoveMap for backbone and sidechain movement
movemap = MoveMap()
movemap.set_bb(True)
movemap.set_chi(True)

# Minimizer setup
min_opts = MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
minimizer = AtomTreeMinimizer()

minimizer.run(mypose, movemap, scorefxn, min_opts)

the_observer = AddPyMOLObserver(mypose)
the_observer.pymol().apply(mypose)

counter_true = 0
counter_false = 0

for i in range(10):
    n_res = mypose.total_residue()
    random_res = numeric.random.uniform()
    residue_index = int(random_res * n_res) + 1

    phi_perturb = numeric.random.gaussian()
    psi_perturb = numeric.random.gaussian()
    orig_phi = mypose.phi(residue_index)
    orig_psi = mypose.psi(residue_index)
    mypose.set_phi(residue_index, orig_phi + phi_perturb)
    mypose.set_psi(residue_index, orig_psi + psi_perturb)

    tf = TaskFactory()
    task = tf.create_task_and_apply_taskoperations(mypose)
    task.restrict_to_repacking()
    pack_rotamers(mypose, scorefxn, task)

    accepted = mc.boltzmann(mypose)

    if accepted == True:
        counter_true += 1
    elif accepted == False:
        counter_false += 1

    if (i+1) % 10 == 0:
        acceptance_rate = counter_true/10
        print(f"Acceptance rate: {acceptance_rate}")
        average_energy = mypose.energies().total_energy() / 10
        print(f"Average energy: {average_energy}")

print(f"Rejections: {counter_false}")
print(f"Acceptances: {counter_true}")

mc.recover_low(mypose)            # overwrite pose with best-so-far
mypose.dump_pdb("out_best.pdb")

print(f"Final score: {scorefxn(mypose)}")

def identify_secondary_structure_spans(ss):
    assert False, "TODO: function that takes a string representing the secondary structure and returns a list describing how many secondary structure elements were found, and who the first and last reisue of each element is"
    elements []
    return elements



