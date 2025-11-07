from pyrosetta import init, pose_from_pdb, get_score_function
from pyrosetta.rosetta import numeric
from pyrosetta.rosetta.protocols.moves import AddPyMOLObserver, MonteCarlo
from pyrosetta.rosetta.core.kinematics import MoveMap, FoldTree
from pyrosetta.rosetta.core.optimization import MinimizerOptions, AtomTreeMinimizer
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack import pack_rotamers
from pyrosetta.rosetta.core.scoring.dssp import Dssp
from pyrosetta.rosetta.core.scoring import ScoreFunction, linear_chainbreak, parse_score_function, attributes_for_parse_score_function_w_description, ScoreType
from pyrosetta.rosetta.core.pose import correctly_add_cutpoint_variants
from pyrosetta.rosetta.protocols.moves import Mover, xsd_type_definition_w_attributes
from pyrosetta.rosetta.utility.tag import XMLSchemaAttribute, XMLSchemaComplexTypeGenerator, XMLSchemaDataType, XMLSchemaCommonType, XMLSchemaType
from pyrosetta.rosetta.std import list_utility_tag_XMLSchemaAttribute_t
from pyrosetta.rosetta import protocols  # for moves.xsd_type_definition_w_attributes
from bootcamp_app import identify_secondary_structure_spans, fold_tree_from_dssp_string


class BootCampMover(Mover):
    
    _clones = list()

    def __init__(self, sfxn: ScoreFunction | None = None, num_iterations: int = 500):
        super().__init__()
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
    
    def get_name(self):
        self.__class__.__name__
        return "BootCampMover"

    def apply(self,pose):
        scorefxn = get_score_function(True)
        current_score = scorefxn(pose)
        print(current_score)
        print("-------------")

        temperature = 1.0
        mc = MonteCarlo(pose, scorefxn, temperature) ##come back and check this syntax

        # Set up MoveMap for backbone and sidechain movement
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)

        # Minimizer setup
        min_opts = MinimizerOptions("lbfgs_armijo_atol", 0.01, True)
        minimizer = AtomTreeMinimizer()

        minimizer.run(pose, movemap, scorefxn, min_opts)

        the_observer = AddPyMOLObserver(pose)
        the_observer.pymol().apply(pose)

        counter_true = 0
        counter_false = 0

        for i in range(10):
            n_res = pose.total_residue()
            random_res = numeric.random.uniform()
            residue_index = int(random_res * n_res) + 1

            phi_perturb = numeric.random.gaussian()
            psi_perturb = numeric.random.gaussian()
            orig_phi = pose.phi(residue_index)
            orig_psi = pose.psi(residue_index)
            pose.set_phi(residue_index, orig_phi + phi_perturb)
            pose.set_psi(residue_index, orig_psi + psi_perturb)

            tf = TaskFactory()
            task = tf.create_task_and_apply_taskoperations(pose)
            task.restrict_to_repacking()
            pack_rotamers(pose, scorefxn, task)

            accepted = mc.boltzmann(pose)

            if accepted == True:
                counter_true += 1
            elif accepted == False:
                counter_false += 1

            if (i+1) % 10 == 0:
                acceptance_rate = counter_true/10
                print(f"Acceptance rate: {acceptance_rate}")
                average_energy = pose.energies().total_energy() / 10
                print(f"Average energy: {average_energy}")

        print(f"Rejections: {counter_false}")
        print(f"Acceptances: {counter_true}")

        mc.recover_low(pose)            # overwrite pose with best-so-far
        pose.dump_pdb("out_best.pdb")

        print(f"Final score: {scorefxn(pose)}")

        ss = Dssp(pose).get_dssp_secstruct()
            # Your helper that finds SSE blocks (expects a string, not a Pose)
        ss1 = identify_secondary_structure_spans(ss)

        # Build FoldTree from DSSP using your helper; support both possible signatures
        try:
            ft = fold_tree_from_dssp_string(pose, ss)  # most likely correct
        except TypeError:
            ft = fold_tree_from_dssp_string(ss) 

        # Build a FoldTree from the DSSP string using your helper
        pose.fold_tree(ft)

        # Add cutpoint variants that correspond to the current FoldTree
        correctly_add_cutpoint_variants(pose)

        # Score with a standard FA scorefunction (or use self._sfxn if you prefer)
        sfxn = get_score_function()
        sfxn.set_weight(ScoreType.linear_chainbreak, 1.0)

        score = sfxn(pose)
        print("Score after fold tree + cutpoints:", score)

    @staticmethod
    def mover_name():
         return "BootCampMover"
    
    def parse_my_tag(self, tag, datamap):
        self.set_scorefunction(parse_score_function(tag, datamap))
        if tag.hasOption("num_iterations"):
            iters = tag.get_option_int("num_iterations", 1)

    @staticmethod
    def provide_xml_schema(xsd):
        attrs = list_utility_tag_XMLSchemaAttribute_t()
        # ---- num_iterations attribute ----
        attrs.append(
            XMLSchemaAttribute.attribute_w_default(
                "num_iterations",
                XMLSchemaType(XMLSchemaCommonType.xsct_positive_integer),
                "Number of iteration",
                "10"
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

    
