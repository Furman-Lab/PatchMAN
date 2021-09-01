#!/vol/ek/Home/alisa/python3.5/bin/python3.5

from pyrosetta import *
from pyrosetta.rosetta import *

import os
import sys

import pyrosetta_utils as utils


MIN_RES_PATCH = 5
NEIGHBORS_DIST = 10
BACKBONE_ATOMS = ['C', 'CA', 'O', 'N']
MAX_STRETCH_LEN = 7
MAX_HELIX_LEN = 11


def define_motifs(pose, pdb_name):

    surf_sel = utils.create_layer_selector()
    selected_res = surf_sel.apply(pose)

    include_focus = True
    neighborhood_selector = utils.create_neighborhood_selector(NEIGHBORS_DIST, include_focus)

    motifs_created = [] # save here res indexes around which patches were created
    motifs_with_chains = [] # residues and chains [('A', '2'), ('A', '3'), ('A', '4')]

    pdbinf = pose.pdb_info()

    all_motifs_list = []
    num = 1
    for i,res in enumerate(selected_res, 1):
        if i <=2:  # skip first 2 res (too small, overlapping patches)
            continue
        if res:
            if i in motifs_created: # skip every second patch (to prevent too overlapping patches)
                continue
            focus = utils.create_index_selector([i])
            neighborhood_selector.set_focus(focus.apply(pose))
            atoms = pyrosetta.rosetta.utility.vector1_std_string()

            atoms.append('CA')

            neighborhood_selector.set_atom_names_for_distance_measure(atoms)
            new_motif_sel = neighborhood_selector.apply(pose)

            motif = list(map(str, core.select.residue_selector.selection_positions(new_motif_sel)))
            refine_motif(motif, pose)
            motif_with_chains = [(pdbinf.chain(int(resn)), str(pdbinf.number(int(resn)))) for resn in motif]

            motifs_created.append(i+1)
            if motif_with_chains not in motifs_with_chains:
                motifs_with_chains.append(motif_with_chains)
                motif_name = '%03d_'%num+pdb_name
                num += 1

                write_to_pdb(motif, motif_name, pose)

        all_motifs_list.append(motif_name)

    return all_motifs_list


def write_to_pdb(motif, motif_name, pose):
    """Dump motif pdbs"""
    motif_pose = Pose()
    slice_res = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for res in motif:
        slice_res.append(int(res))
    core.pose.pdbslice(motif_pose, pose, slice_res)
    motif_pose.dump_pdb(motif_name + '.pdb')


def refine_motif(motif, pose):
    """Stretches are not longer than MAX_LENGTH, no single residues"""
    stretch_len = 1
    res_to_remove = []
    pose.display_secstruct()
    s = pose.secstruct()
    for idx, resid in enumerate(motif):
        if idx == len(motif) - 1:  # The END
            break
        elif int(resid) + 1 == int(motif[idx + 1]):
            stretch_len += 1
        elif stretch_len > 1:
            if is_helix(motif[idx + 1 - stretch_len:idx + 1], s):
                max_length = MAX_HELIX_LEN
            else:
                max_length = MAX_STRETCH_LEN
            if stretch_len > max_length:
                res_to_remove.extend(motif[idx - (stretch_len - max_length) + 1:idx + 1]) # cut long sequences
            stretch_len = 1
            continue
        else:
            res_to_remove.append(resid)
    if stretch_len <= 1:
        res_to_remove.append(resid)
    elif stretch_len > MAX_STRETCH_LEN and stretch_len == len(motif): # cut long sequences
        if is_helix(motif[idx + 1 - stretch_len:idx + 1], s):
            max_length = MAX_HELIX_LEN
        else:
            max_length = MAX_STRETCH_LEN
        res_to_remove = motif[max_length:]
    for r in res_to_remove:
        motif.remove(r)


def is_helix(stretch, s):
    if s[int(stretch[0]):int(stretch[-1]) + 1].count('H') >= 3:
        return True
    else:
        return False


def main():

    inpdb = sys.argv[1] # CLEAN_PDB
    toolbox.cleaning.cleanATOM(inpdb)
    prot_name = os.path.splitext(os.path.basename(inpdb))[0]
    pose = pose_from_pdb(prot_name+'.clean.pdb')
    motifs = define_motifs(pose,prot_name)
    print("The surface was split into " + str(len(motifs)) + " patches")


if __name__ == "__main__":

    pyrosetta.init()

    main()
