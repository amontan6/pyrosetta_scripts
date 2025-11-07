#PATH/TO/ROSETTA/source/src/core/simple_metrics/per_residue_metric
from pyrosetta import *
pyrosetta.init()
from pyrosetta.rosetta.utility.tag import XMLSchemaAttribute, XMLSchemaType
from pyrosetta.rosetta.utility.tag import xs_string
from pyrosetta.rosetta.core.simple_metrics import PerResidueRealMetric, xsd_simple_metric_type_definition_w_attributes
from pyrosetta.rosetta.std import list_utility_tag_XMLSchemaAttribute_t
from pyrosetta.rosetta.core.pose import PDBInfo

#pose = pose_from_pdb("AdenosineRec_AF.pdb")

class PerResidueBfactorBootCampMetric(PerResidueRealMetric):
    def _init_(self, atom_type_):
        super().__init__()
        self._atom_type_ = atom_type_
    
    _clones = list()

    def calculate(self):
        pdb_info = pose.pdb_info()
        bfactors = {}
        for i in range(1, pose.total_residue() + 1):
            atom_b = [pdb_info.bfactor(i, j) for j in range(1, pose.residue(i).natoms() + 1)]
            bfactors[i] = sum(atom_b) / len(atom_b)
        return bfactors

    def name(self) -> str:
        self.__class__.__name__
        return "PerResidueBfactorBootCampMetric"

    def metric(self, B_factormetric) -> str:
        return B_factormetric ###???? come back to this!

    def parse_my_tag(self,tag, datamap):
        if tag.hasOption("atom_type_"):
            iters = tag.get_option_string("atom_type_")

    def clone(self):
        copy = PerResidueBfactorBootCampMetric()
        copy.atom_type_ = self.atom_type_
        PerResidueBfactorBootCampMetric._clones.append(copy)
        return copy
    
    @classmethod
    def provide_xml_schema(cls, xsd):
        attrlist = list_utility_tag_XMLSchemaAttribute_t()

        attrlist.append(XMLSchemaAttribute.required_attribute(
                "atom_type_",
                XMLSchemaType(xs_string),
                "FILL OUT THIS DESCRIPTION"))

        description = '''
            FILL OUT THIS DESCTIPTION 
    '''
        xsd_simple_metric_type_definition_w_attributes(
                xsd,
                cls.name(),
                description, attrlist)