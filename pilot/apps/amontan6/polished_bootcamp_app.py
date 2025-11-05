#!/usr/bin/env python3
import argparse
from dataclasses import dataclass

# -------- PyRosetta imports (safe to import; no init yet) --------
from pyrosetta import init, pose_from_pdb, get_score_function
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.kinematics import MoveMap, FoldTree
from pyrosetta.rosetta.core.optimization import MinimizerOptions, AtomTreeMinimizer
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import linear_chainbreak
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants


# ==============================
# Helpers (import-safe)
# ==============================
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

# ==============================
# FoldTree builders (import-safe)
# ==============================
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
    sse  = identify_secondary_structure_spans(ss)
    if not sse:
        raise ValueError("No SSEs (H/E) found in DSSP string.")
    loops = _loops_from_sse(sse)

    sse_m  = [_mid(a, b) for (a, b) in sse]
    loop_m = [_mid(a, b) for (a, b) in loops]
    root   = sse_m[0]

    ft = FoldTree()
    jump_id = 1

    # Root SSE peptides (back to 1, forward to end of SSE1)
    a0, b0 = sse[0]; m0 = sse_m[0]
    if m0 > 1:
        ft.add_edge(m0, 1, -1)     # (root → 1)
    if m0 < b0:
        ft.add_edge(m0, b0, -1)    # (root → end of SSE1)

    # Interleave rows: loop_k then SSE_{k+1}
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
            forward_to = N if (k + 1) == (len(sse) - 1) else sb
            if sm < forward_to: ft.add_edge(sm, forward_to, -1)

    # Optional: root the tree (edge order unchanged)
    try:
        ft.reorder(root)
    except Exception:
        pass

    return ft

def fold_tree_from_ss(pose) -> FoldTree:
    """Compute DSSP on the pose and return a FoldTree built from the DSSP string."""
    ss = Dssp(pose).get_dssp_secstruct()
    return fold_tree_from_dssp_string(pose, ss)

# ==============================
# CLI + main (run-time only)
# ==============================
def build_parser():
    p = argparse.ArgumentParser(description="PyRosetta FoldTree demo")
    p.add_argument("filename", help="Path to the input PDB file")
    p.add_argument("--trials", type=int, default=10, help="MC trials (default: 10)")
    p.add_argument("--temperature", type=float, default=1.0, help="MC temperature")
    p.add_argument("--observe", action="store_true", help="Attach PyMOL observer")
    return p

def main(argv=None):
    args = build_parser().parse_args(argv)

    # Initialize Rosetta
    init(extra_options="-ignore_unrecognized_res")

    # Load pose
    pose = pose_from_pdb(args.filename)
    print(f"Loaded pose with {pose.total_residue()} residues from: {args.filename}")

    # Build & apply FoldTree from DSSP
    ft = fold_tree_from_ss(pose)
    print(ft)
    pose.fold_tree(ft)

    # Enable chainbreak handling (must be AFTER setting the FoldTree)
    scorefxn = get_score_function(True)
    scorefxn.set_weight(linear_chainbreak, 1.0)
    correctly_add_cutpoint_variants(pose)

    # Score before moves
    print("Initial score:", scorefxn(pose))
    print("-------------")

    # Optional PyMOL observer
    if args.observe:
        obs = AddPyMOLObserver(pose)
        obs.pymol().apply(pose)

    # MoveMap & minimizer
    movemap = MoveMap(); movemap.set_bb(True); movemap.set_chi(True)
    min_opts = MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
    AtomTreeMinimizer().run(pose, movemap, scorefxn, min_opts)

    # Monte Carlo loop
    mc = MonteCarlo(pose, scorefxn, args.temperature)
    accepts = rejects = 0
    for i in range(args.trials):
        n_res = pose.total_residue()
        residue_index = int(numeric.random.uniform() * n_res) + 1

        #debugging to see energies
        E_before = scorefxn(pose)

        # Random small torsion tweaks
        phi_d  = numeric.random.gaussian()
        psi_d  = numeric.random.gaussian()
        pose.set_phi(residue_index, pose.phi(residue_index) + phi_d)
        pose.set_psi(residue_index, pose.psi(residue_index) + psi_d)

        # Repack sidechains
        tf = TaskFactory()
        task = tf.create_task_and_apply_taskoperations(pose)
        task.restrict_to_repacking()
        pack_rotamers(pose, scorefxn, task)

        # Energy after move
        E_after = scorefxn(pose)
        dE = E_after - E_before

        if mc.boltzmann(pose):
            accepts += 1
        else:
            rejects += 1

        print(f"Trial {i+1:3d}: res {residue_index:3d}, ΔE = {dE:8.3f}, "
          f"E_before = {E_before:8.3f}, E_after = {E_after:8.3f}, "
          f"accepted = {accepts}")

    print(f"Rejections: {rejects}")
    print(f"Acceptances: {accepts}")
    print(f"Final acceptance rate: {accepts / args.trials}")
    mc.recover_low(pose)
    pose.dump_pdb("out_best.pdb")
    print(f"Final score: {scorefxn(pose)}")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())

"""One of the issues with the fact that that the acceptances are no low, or 0, 
is because I am not doing something to fix the breaks, I'm just breaking the chain?? \(-_-)/ """
