from pyrosetta import init, pose_from_pdb, get_score_function
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.kinematics import MoveMap, FoldTree
from pyrosetta.rosetta.core.optimization import MinimizerOptions, AtomTreeMinimizer
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import ScoreFunction, linear_chainbreak, parse_score_function, attributes_for_parse_score_function_w_description
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants
from pyrosetta.rosetta.protocols.moves import Mover

class BootCampMover(Mover):
    
    _clones = list()

    def __init__(self, sfxn: ScoreFunction | None = None, num_iterations: int = 500):
        super().__init__(self)
        self._sfxn = sfxn if sfxn is not None else get_score_function()
        self._num_iterations = num_iterations

    def clone(self):
        copy = BootCampMover(self._sfxn, self._num_iterations)
        BootCampMover._clones.append(copy)
        return copy

    def fresh_instance(self):
        return BootCampMover()

    def set_scorefunction(self, sfxn: ScoreFunction) -> None:
        self._sfxn = sfxn

    def scorefunction(self) -> ScoreFunction:
        return self._sfxn

    def set_num_iterations(self, n: int) -> None:
        self._num_iterations = int(n)

    def num_iterations(self) -> int:
        return self._num_iterations

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

    def mover_name(self):
         return "BootCampMover"
    
    def parse_my_tag(self, tag, datamap):
        if tag.hasOption("num_iterations"):
            iters = tag.get_option_int("num_iterations", 1)
        parse_score_function(self._sfxn)

    @staticmethod
    def provide_xml_schema(xsd):
        attrs = list_utility_tag_XMLSchemaAttribute_t()
        # ---- num_iterations attribute ----
        attrs.append(
            XMLSchemaAttribute.attribute_w_default(
                "num_iterations",
                XMLSchemaCommonType.xsct_positive_integer,   # required type
                "10",                                       # default
                "Number of Monte Carlo refinement iterations."
            )
        )
        # ---- ScoreFunction attribute helper ----
        attributes_for_parse_score_function_w_description(
            attrs,
            "ScoreFunction to use for sampling"
        )
        # Register this mover’s schema
        xsd_type_definition_w_attributes(
            xsd,
            BootCampMover.mover_name(),
            "BootCampMover constructs a FoldTree from DSSP and performs refinement using Monte Carlo.",
            attrs
        )

    
