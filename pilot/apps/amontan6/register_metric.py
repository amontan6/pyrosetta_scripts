from pyrosetta.rosetta.core.simple_metrics import SimpleMetricCreator, SimpleMetricFactory
from Per_Residue_Bfactor_Metric import PerResidueBfactorBootCampMetric

class PerResidueBfactorBootCampMetricCreator(SimpleMetricCreator):
    instances_ = list()
    _py_mover_creators_ = []

    def __init__(self):
        SimpleMetricCreator.__init__(self)

    def create_simple_metric(self):
        B_factormetric = PerResidueBfactorBootCampMetric()
        self.instances_.append(B_factormetric)
        return B_factormetric

    def keyname(self):
        return PerResidueBfactorBootCampMetric.name(self)

    def provide_xml_schema(self, xsd):
        print("creator provide_xml_schema is called")
        PerResidueBfactorBootCampMetric.provide_xml_schema(xsd)

    @staticmethod
    def register():
        factory = SimpleMetricFactory.get_instance()
        creator = PerResidueBfactorBootCampMetricCreator()
        factory.factory_register(creator)

        PerResidueBfactorBootCampMetricCreator._py_mover_creators_.append(creator)

