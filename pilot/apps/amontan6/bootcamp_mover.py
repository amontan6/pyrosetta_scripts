import pyrosetta.rosetta.protocols.moves import Mover
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

class BootCampMover(Mover):

    def __init__(self):
        super().__init__(self)

    def apply(self,pose):
            ss = Dssp(pose).get_dssp_secstruct()
            ft = fold_tree_from_dssp_string(ss)
            pose.fold_tree(ft)
            correctly_add_cutpoint_variants(pose)
            sfxn = get_score_function()
            sfxn.set_weight(linear_chainbreak, 1.0)
            score = sfxn(pose)
            print("Score after fold tree + cutpoints:", score)

    def get_name(self):
        self.__class__.__name__
