#!/vol/ek/Home/alisa/python3.5/bin/python3

from pyrosetta import *
from pyrosetta.rosetta import *
from rosetta_protocols import prepack
import sys

init('-ex1 -ex2aro -use_input_sc')

with open('all_pdbs', 'r') as fh:
    all_pdbs = list(set([l.strip() for l in fh.readlines()]))

for p in all_pdbs:
    pdb_pose = toolbox.pose_from_rcsb(p)

rec = sys.argv[1]
rec_name = rec[:-4]
pose = pose_from_pdb(sys.argv[1])
prepack(pose)

pose.dump_pdb(rec_name + '.ppk.pdb')
