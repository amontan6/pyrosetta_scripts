import sys
import argparse
from pyrosetta import *
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.kinematics import MoveMap, FoldTree
from pyrosetta.rosetta.core.optimization import MinimizerOptions, AtomTreeMinimizer
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from dataclasses import dataclass

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

#dssp = Dssp(mypose)
#dssp_calc = dssp.get_dssp_secstruct()
#print(dssp_calc)

#ss = ""
#print(identify_secondary_structure_spans(dssp_calc))

def identify_secondary_structure_spans(ss: str):
    """Contiguous runs of the SAME H or E (1-based inclusive spans)."""
    blocks, start, current = [], None, None
    for i, ch in enumerate(ss, start=1):
        if ch in ('E', 'H'):
            if current is None:
                current, start = ch, i
            elif ch != current:
                blocks.append((start, i - 1))
                current, start = ch, i
        else:
            if current is not None:
                blocks.append((start, i - 1))
                current, start = None, None
    if current is not None:
        blocks.append((start, len(ss)))
    return blocks

def _loops_from_sse(sse_blocks):
    """Inter-SSE loops = non-H/E gaps between adjacent SSEs (no terminal loops)."""
    loops = []
    for i in range(len(sse_blocks) - 1):
        a_end = sse_blocks[i][1]
        b_beg = sse_blocks[i+1][0]
        if b_beg - a_end > 1:
            loops.append((a_end + 1, b_beg - 1))
    return loops

def _mid(a: int, b: int) -> int:
    """1-based floor midpoint."""
    return (a + b) // 2


# ---------- core: build FoldTree from DSSP string (doc order) ----------
def fold_tree_from_dssp_string(pose, ss: str) -> FoldTree:
    """
    Build a FoldTree in the exact document order:
      root back, root forward,
      then for k = 0..:  Jump(root→loop_k), loop back, loop fwd,
                         Jump(root→sse_{k+1}), sse back, sse fwd
    Special cases:
      - root back goes to residue 1 (not just SSE start)
      - last SSE forward goes to residue N (pose length)
    """
    if not ss:
        raise ValueError("Empty secondary-structure string.")
    if len(ss) != pose.total_residue():
        raise ValueError(f"SS string length {len(ss)} != pose length {pose.total_residue()}")

    N    = pose.total_residue()
    sse  = identify_secondary_structure_spans(ss)   # H/E blocks (same-letter runs)
    if not sse:
        raise ValueError("No SSEs (H/E) found in DSSP string.")
    loops = _loops_from_sse(sse)                    # non-H/E gaps between SSEs

    sse_m  = [_mid(a, b) for (a, b) in sse]
    loop_m = [_mid(a, b) for (a, b) in loops]
    root   = sse_m[0]

    ft = FoldTree()
    jump_id = 1

    # --- Root SSE peptides (back to 1, forward to end of SSE1) ---
    a0, b0 = sse[0]; m0 = sse_m[0]
    if m0 > 1:
        ft.add_edge(m0, 1, -1)     # e.g., (7 → 1)
    if m0 < b0:
        ft.add_edge(m0, b0, -1)    # e.g., (7 → 10)

    # --- Interleave rows: loop_k then SSE_{k+1} ---
    rows = max(len(loops), len(sse) - 1)
    for k in range(rows):
        # loop_k: Jump, back, forward (if exists)
        if k < len(loops):
            lm = loop_m[k]; la, lb = loops[k]
            ft.add_edge(root, lm, jump_id); jump_id += 1
            if lm > la: ft.add_edge(lm, la, -1)
            if lm < lb: ft.add_edge(lm, lb, -1)

        # sse_{k+1}: Jump, back, forward (if exists)
        if (k + 1) < len(sse):
            sm = sse_m[k+1]; sa, sb = sse[k+1]
            ft.add_edge(root, sm, jump_id); jump_id += 1
            if sm > sa: ft.add_edge(sm, sa, -1)
            # forward: last SSE goes all the way to N
            forward_to = N if (k + 1) == (len(sse) - 1) else sb
            if sm < forward_to: ft.add_edge(sm, forward_to, -1)

    #root the tree at 'root' (order of edges is unchanged)
    try:
        ft.reorder(root)
    except Exception:
        pass

    return ft


# ---------- convenience: Pose → DSSP → FoldTree ----------
def fold_tree_from_ss(pose) -> FoldTree:
    """
    Compute DSSP on the pose and return a FoldTree built from the DSSP string.
    """
    ss = Dssp(pose).get_dssp_secstruct()
    return fold_tree_from_dssp_string(pose, ss)

ft = fold_tree_from_ss(mypose)
print(ft)