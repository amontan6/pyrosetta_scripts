@dataclass
class _Built:
    ft: FoldTree
    loops: List[Loop]
    loop_for_residue: List[int]  # 1..N, values in 0..len(loops)

class FoldTreeFromSS:

   def __init__(self, pose: core.pose.Pose, loop_left: int = 2, loop_right: int = 3):
   
   def fold_tree(self) -> FoldTree:
       assert False, "TODO"


   def loop(self, index: int) -> Loop:
       assert False, "TODO"

   def loop_for_residue(self, seqpos: int) -> int:
       assert False, "TODO"