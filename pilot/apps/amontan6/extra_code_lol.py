"""
def fold_tree_from_dssp_string(pose, ss_string, loop_cut_min=2):
    '''Outputs Fold Tree of Inputted Pose'''
    # Total number of residues that the fold tree must span
    N = pose.total_residue()
    # Intialize FoldTree with length N
    ft = FoldTree()
    
    # Identify long loop regions which would make a cycle so they need to be cut
    loop_runs = [] # Keeps track of long loop start and stop indicies
    run_start = None
    for i,s in enumerate(ss_string, start=1):
        # Checks to see if certain parts of the string are in a flexible loop region
        is_loop = (s not in {'H', 'E'})
        if is_loop:
            if run_start is None:
                # If it is we are going to record the start index of the loop
                run_start = i
        else:
            # Checks to see what the stop index is
            if run_start is not None:
                # Check to see if the loop is long (basically if its more than one L)
                if (i - run_start) >= loop_cut_min:
                    # If it is a long loop, we add the start and stop index to loop_runs
                    loop_runs.append((run_start, i-1))
                # Set run_start back to None to restart start and stop block
                run_start = None
    # Records final loop near C-terminus
    if run_start is not None and (N+1 - run_start) >= loop_cut_min:
        loop_runs.append((run_start,N))
    # Walk though each loop to construct full fold tree
    jump_id = 1
    last_anchor = 1 # Tracks where last segment started
    # H H H(loop_start) L L L(loop_end)H H H
    for (loop_start, loop_end) in loop_runs:
        # Peptide edge end point
        seg_end = loop_start - 1
        if last_anchor < seg_end:
            # Record peptide start and stop
            ft.add_edge(last_anchor, seg_end,-1)
        up = seg_end
        down = loop_end + 1
        if 1 <= up <= N and 1 <= down <= N and up < down:
            # Add jump across loop region
            ft.add_edge(up,down, jump_id)
            jump_id +=1
        last_anchor = loop_end+1
    # Adds final C-terminus peptide edge
    if last_anchor < N:
        ft.add_edge(last_anchor,N,-1)
    print(ft)
    pose.fold_tree(ft)
    return ft
"""
"""
def fold_tree_from_dssp_string(pose, ss_string, loop_cut_min=2):
    N = pose.total_residue()
    if len(ss_string) != N:
        raise ValueError(f"SS string length {len(ss_string)} != pose length {N}")

    # Start from an empty tree (no default 1–N edge)
    ft = FoldTree()

    # --- find loop runs: anything not H/E is loop/coil ---
    loop_runs = []
    run_start = None
    for i, s in enumerate(ss_string, start=1):
            if s not in {'H', 'E'}:
                if run_start is None:
                    run_start = i
            else:
                if run_start is not None:
                    if (i - run_start) >= loop_cut_min:
                        loop_runs.append((run_start, i - 1))
                    run_start = None
    if run_start is not None and (N + 1 - run_start) >= loop_cut_min:
        loop_runs.append((run_start, N))

    # --- build edges/jumps with a special case for a C-terminal loop ---
    jump_id = 1
    last_anchor = 1

    for (loop_start, loop_end) in loop_runs:
        seg_end = loop_start - 1

        # C-terminal loop: do NOT cut; connect through to N and finish.
        if loop_end == N:
            if last_anchor < N:
                ft.add_edge(last_anchor, N, -1)   # peptide all the way to C-ter
            last_anchor = N + 1
            # We break here because nothing exists beyond a C-terminal loop.
            break

    # Internal loop: peptide before it, then jump over it
        if last_anchor < seg_end:
            ft.add_edge(last_anchor, seg_end, -1)      # peptide
        up, down = seg_end, loop_end + 1
        if 1 <= up <= N and 1 <= down <= N and up < down:
            ft.add_edge(up, down, jump_id)             # jump
            jump_id += 1
        last_anchor = loop_end + 1

    # If we never hit a terminal loop and there’s trailing structured segment, add it
    if last_anchor < N:
        ft.add_edge(last_anchor, N, -1)

    pose.fold_tree(ft)
    return ft

def fold_tree_from_ss(pose, loop_cut_min=2):
    dssp = Dssp(pose)
    ss_string = dssp.get_dssp_secstruct()
    print(ss_string)
    ft = fold_tree_from_dssp_string(pose, ss_string, loop_cut_min)
    return ft

print(fold_tree_from_ss(mypose))

#def fold_tree_from_ss(pose):
#   ss = Dssp(pose).get_dssp_secstruct()
#   return fold_tree_from_dssp_string(pose, ss)
"""
"""
# ---- helpers ----
def _mid(a: int, b: int) -> int:
    return (a + b) // 2  # 1-based floor midpoint

def _loops_from_sse(sse_blocks):
    loops = []
    for i in range(len(sse_blocks) - 1):
        a_end = sse_blocks[i][1]
        b_beg = sse_blocks[i+1][0]
        if b_beg - a_end > 1:
            loops.append((a_end + 1, b_beg - 1))
    return loops

def pretty_fold_tree(edges):
    parts = ["FOLD_TREE"]
    for u, v, lab in edges:
        parts.append(f" EDGE {u} {v} {lab}")
    return " ".join(parts)

# ---- STRING-ONLY: builds the edges in the exact doc order you described ----
def fold_tree_from_dssp_string(ss: str):

    Returns a list of (u, v, label) edges in the *document order*:
      root back, root forward,
      [for each row k] Jump(root→loop_k), loop back, loop fwd,
                       Jump(root→sse_{k+1}), sse back, sse fwd
    label = -1 for peptide, >0 for jump.

    if not ss:
        raise ValueError("Empty secondary-structure string.")

    sse = identify_secondary_structure_spans(ss)   # your parser for H/E blocks
    if not sse:
        raise ValueError("No SSEs (H/E) found.")
    loops = _loops_from_sse(sse)                   # gaps between SSEs (no terminal loops)

    sse_m  = [_mid(a, b) for (a, b) in sse]
    loop_m = [_mid(a, b) for (a, b) in loops]
    root   = sse_m[0]
    N      = len(ss)

    edges = []
    jump_id = 1

    # Root SSE peptides: back to residue 1, then forward to end of SSE1
    a0, b0 = sse[0]; m0 = sse_m[0]
    if m0 > 1:
        edges.append((m0, 1, -1))       # e.g., (7 → 1)
    if m0 < b0:
        edges.append((m0, b0, -1))      # e.g., (7 → 10)

    # Interleave rows: loop_k then SSE_{k+1}
    rows = max(len(loops), len(sse) - 1)
    for k in range(rows):
        # loop_k block
        if k < len(loops):
            lm = loop_m[k]; la, lb = loops[k]
            edges.append((root, lm, jump_id)); jump_id += 1      # Jump(root→loop_k)
            if lm > la: edges.append((lm, la, -1))               # loop back
            if lm < lb: edges.append((lm, lb, -1))               # loop forward

        # SSE_{k+1} block
        if (k + 1) < len(sse):
            sm = sse_m[k+1]; sa, sb = sse[k+1]
            edges.append((root, sm, jump_id)); jump_id += 1      # Jump(root→SSE_{k+1})
            if sm > sa: edges.append((sm, sa, -1))               # SSE back
            # 🔧 Forward: if this is the LAST SSE, go to N; else go to its own end
            forward_to = N if (k + 1) == (len(sse) - 1) else sb
            if sm < forward_to:
                edges.append((sm, forward_to, -1))               # SSE forward (last goes to N)

    return edges


example = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     "
ft = fold_tree_from_dssp_string(example)
print(ft)
"""

# ---------- helpers ----------
