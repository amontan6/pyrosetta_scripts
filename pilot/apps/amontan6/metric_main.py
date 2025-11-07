from pyrosetta import *
from pyrosetta.rosetta.protocols.moves import MoverFactory
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from Per_Residue_Bfactor_Metric import PerResidueBfactorBootCampMetric
from register_metric import PerResidueBfactorBootCampMetricCreator

EMBEDDED_XML = """
<ROSETTASCRIPTS>
	<SCOREFXNS>
	</SCOREFXNS>
	</RESIDUE_SELECTORS>
	<MOVE_MAP_FACTORIES>
	</MoveMapFactory>
	</MOVE_MAP_FACTORIES>
	<SIMPLE_METRICS>
        <PerResidueBfactorBootCampMetric name="bfactormetric" residue_selector="ALL"/>
	</SIMPLE_METRICS>
	<MOVERS>
	</MOVERS>
	<PROTOCOLS>
        <Add mover_name="bfactormetric"/>
	</PROTOCOLS>
</ROSETTASCRIPTS>
"""
def main():
    PerResidueBfactorBootCampMetricCreator.register()
    init()
    pose = pose_from_pdb("AdenosineRec_AF.pdb ")
    
    xmlobj = XmlObjects.create_from_string(EMBEDDED_XML)
    protocol = xmlobj.get_mover("ParsedProtocol")
    protocol.apply(pose)

    pose.dump_pdb("Bfactor_output.pdb")
    print(" Finished! Output written to Bfactor_output.pdb")

if __name__ == "__main__":
    main()