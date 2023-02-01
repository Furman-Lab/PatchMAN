

from pyrosetta import *
from pyrosetta.rosetta import *

import os
import sys

import pyrosetta_utils as utils
import argparse


MIN_RES_PATCH = 5
NEIGHBORS_DIST = 10
BACKBONE_ATOMS = ['C', 'CA', 'O', 'N']
MAX_STRETCH_LEN = 7
MAX_HELIX_LEN = 11
def define_motifs(pose, pdb_name, mask_pose=None, focus_mode=False):

    surf_sel = utils.create_layer_selector()
    selected_res = surf_sel.apply(pose)
    include_focus = True
    neighborhood_selector = utils.create_neighborhood_selector(NEIGHBORS_DIST, include_focus)

    motifs_created = [] # save here res indexes around which patches were created
    motifs_with_chains = [] # residues and chains [('A', '2'), ('A', '3'), ('A', '4')]

    pdbinf = pose.pdb_info()
    if mask_pose is not None:
        mask_pose_pdbinf = mask_pose.pdb_info()

    all_motifs_list = []
    num = 1
    focus_res_motif = []
    for i, res in enumerate(selected_res, 1):
        if i <= 2:  # skip first 2 res (too small, overlapping patches)
            continue
        if res:
            if i in motifs_created: # skip every second patch (to prevent too overlapping patches)
                continue

            # if focus mode is on, do not use residues as patch center that are not in the focus pose
            focus_res_chain = pdbinf.chain(int(i))
            focus_res_number = pdbinf.number(int(i))
            if focus_mode:
                if mask_pose_pdbinf.pdb2pose(focus_res_chain, focus_res_number) == 0:
                    continue
            focus = utils.create_index_selector([i])
            neighborhood_selector.set_focus(focus.apply(pose))
            atoms = pyrosetta.rosetta.utility.vector1_std_string()
            atoms.append('CA')
            neighborhood_selector.set_atom_names_for_distance_measure(atoms)
            new_motif_sel = neighborhood_selector.apply(pose)

            motif = list(map(str, core.select.residue_selector.selection_positions(new_motif_sel)))
            refine_motif(motif, pose)

            motif_not_masked = True
            number_of_masked_residues = 0
            if mask_pose is not None:
                motif_with_chains = []
                motif_length = len(motif)
                
                required_unmasked_length = motif_length * RATIO_ALLOWED_MASKED_RESIDUES
                for resn in motif: # checking for masked residues in the motif
                    resi_in_mask = mask_pose_pdbinf.pdb2pose(motif_res_chain, motif_res_number)
                    motif_res_chain = pdbinf.chain(int(resn))
                    motif_res_number = pdbinf.number(int(resn))

                    if resi_in_mask == 0: # if pdb chain and residue are not in mask, add to motif
                        motif_with_chains.append((motif_res_chain, str(motif_res_number)))
                    else:
                        number_of_masked_residues += 1
                    
                    # break the loop if too many residues are masked
                    if number_of_masked_residues > required_unmasked_length:
                        motif_not_masked = False
                        break
                        
                ratio_of_masked = round(number_of_masked_residues/motif_length, 2)
                print("{masked}/{motiflen} =  {ratio} residues in the motif are masked. Motif kept: {kept}".format(motiflen=motif_length,
                                                    masked=number_of_masked_residues, ratio=ratio_of_masked, kept=motif_not_masked))
            else:
                motif_with_chains = [(pdbinf.chain(int(resn)), str(pdbinf.number(int(resn)))) for resn in motif]

            if motif_not_masked: # add checking if motif contained masked residues
                motifs_created.append(i+1)
                if motif_with_chains not in motifs_with_chains:
                    motifs_with_chains.append(motif_with_chains)
                    motif_name = '%03d_'%num+pdb_name
                    num += 1

                    write_to_pdb(motif, motif_name, pose)
                    all_motifs_list.append(motif_name) # is this the right level of indentation?
                    focus_res_motif.append(['%03d' % num, str(pdbinf.number(int(i))), ','.join(motif)])
        
    with open('focus_res_motif.txt', 'w') as f:
        for l in focus_res_motif:
            f.write(' '.join(l) + '\n')

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
    #pose.display_secstruct()
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
    
def load_and_clean_pdb(pdb_file, return_name=False):
    toolbox.cleaning.cleanATOM(pdb_file)
    prot_name = os.path.splitext(os.path.basename(pdb_file))[0]
    pose = pose_from_pdb(prot_name+'.clean.pdb')
    
    if return_name:
        return pose, prot_name
    else:
        return pose

def main():
    # set up 3 arguments with argparse: -i, -f and -m
    parser = argparse.ArgumentParser(description='Create motifs from a protein structure')
    parser.add_argument('-i', '--input', help='Input PDB file', required=True)
    parser.add_argument('-f', '--focus', help='Focus mask PDB file, with the residues that should form the interface', required=False, default=None)
    parser.add_argument('-m', '--mask', help='Mask PDB file, with residues that should not be on the interface', required=False, default=None)
    
    args = parser.parse_args()
    
    if args.mask is not None and args.focus is not None:
        print('ERROR: Please provide either a mask or a focus file, not both!')
        exit(1)
    
    # load receptor pose
    pose, prot_name = load_and_clean_pdb(args.input, return_name=True)

    if args.mask is not None: # mask PDB file, containing residues NOT to use
        print('reading masking residues from: ' + args.mask )
        mask_pose = load_and_clean_pdb(args.mask)
        motifs = define_motifs(pose, prot_name, mask_pose)
    elif args.focus is not None: # focus PDB file, containing residues to be included in the interface
        print('reading focus residues from: ' + args.focus )
        focus_pose = load_and_clean_pdb(args.focus)
        motifs = define_motifs(pose, prot_name, focus_pose, focus_mode=True)
    else:
        motifs = define_motifs(pose, prot_name)

    print("The surface was split into " + str(len(motifs)) + " patches")


if __name__ == "__main__":

    pyrosetta.init('-mute warning -mute info -mute debug -mute trace -mute core.select.residue_selector.NeighborhoodResidueSelector -mute core')

    main()


RATIO_ALLOWED_MASKED_RESIDUES = 0.8 # added for masking
