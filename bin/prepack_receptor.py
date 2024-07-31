

import sys

from pyrosetta import *
from pyrosetta.rosetta import *

from rosetta_protocols import prepack


def main():

    rec = sys.argv[1]

    init('-ex1 -ex2aro -use_input_sc -unboundrot %s'%rec)  # for pyrosetta

    rec_pose = pose_from_pdb(rec)

    ppk_pose = prepack(rec_pose)

    ppk_pose.dump_pdb(rec[:-4] + '.ppk.pdb')


if __name__ == "__main__":

    main()