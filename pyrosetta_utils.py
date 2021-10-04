#!/vol/ek/Home/alisa/python3.5/bin/python3.5

from pyrosetta import *
from pyrosetta.rosetta import *
import numpy as np


def create_layer_selector():
    lay_sel = rosetta.core.select.residue_selector.LayerSelector()
    lay_sel.set_layers(0, 0, 1)
    lay_sel.set_use_sc_neighbors(False)
    lay_sel.set_ball_radius(1.35)
    return lay_sel


def create_index_selector(res_nums):
    idx_selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
    for res in res_nums:
        idx_selector.append_index(res)
    return idx_selector


def create_neighborhood_selector(cutoff, include_focus):

    neighborhood_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    neighborhood_selector.set_include_focus_in_subset(include_focus)
    neighborhood_selector.set_distance(cutoff)

    return neighborhood_selector


def two_atoms_distance(complex_pose, res1, atom1, res2, atom2):
    coord1 = complex_pose.residue(res1).atom(atom1).xyz()
    coord2 = complex_pose.residue(res2).atom(atom2).xyz()
    x1, y1, z1 = coord1
    x2, y2, z2 = coord2
    distance = np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)
    return distance
