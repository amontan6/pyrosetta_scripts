# Embedded RosettaScripts XML:
from pyrosetta import init, pose_from_pdb
from pyrosetta.rosetta.protocols.rosetta_scripts import XmlObjects
from register_mover import register

# 1) Register before init so factory knows your mover at startup.
register()
init()

# 2) Read input pose
pose = pose_from_pdb("input.pdb")

EMBEDDED_XML = """
<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="sfxn" weights="ref2015"/>
  </SCOREFXNS>
  <MOVERS>
  </MOVERS>
    <BootCampMover, name="bcm" num_iterations="" scorefxn="sfxn"/>
  <PROTOCOLS>
    <Add filter_name="bcm"/>
  </PROTOCOLS>
</ROSETTASCRIPTS>
"""

xmlobj = XmlObjects.create_from_string(EMBEDDED_XML)
protocol = xmlobj.get_mover("ParsedProtocol")
protocol.apply(pose)

pose.dump_pdb("bootcamp_out.pdb")