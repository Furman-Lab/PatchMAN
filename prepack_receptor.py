#!/vol/ek/Home/alisa/python3.5/bin/python3

import sys

from pyrosetta import *
from pyrosetta.rosetta import *

from rosetta_protocols import prepack


def main():
    rec_pose = pose_from_pdb(sys.argv[1])

    ppk_pose = prepack(rec_pose)

    ppk_pose.dump_pdb(sys.argv[1][:-4] + '.ppk.pdb')


if __name__ == "__main__":

    init('-ex1 -ex2aro -ex3 -ex4 -use_input_sc -unboundrot 1Y8T_A.pdb')  # for pyrosetta

    main()