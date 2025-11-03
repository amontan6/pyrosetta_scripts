import sys
import argparse
from pyrosetta import *

init(extra_options="-ignore_unrecognized_res")

parser = argparse.ArgumentParser(description="PyRosetta file loader example")

parser.add_argument("filename", help="Path to the input PDB file")
args = parser.parse_args()

print(f"Filename: {args.filename}")

mypose = pose_from_pdb(args.filename)

print(f"Loaded pose with {mypose.total_residue()} residues from: {args.filename}")