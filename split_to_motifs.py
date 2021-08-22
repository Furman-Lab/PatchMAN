#!/vol/ek/Home/alisa/python3.5/bin/python3.5

# import common_func as utils
# import geohash_create_table as geohash
# import rosetta_protocols as rp
# import hashing as hsh

from pyrosetta import *
from pyrosetta.rosetta import *

import os
import argparse
import pickle
import pandas as pd
import multiprocessing as mp
from typing import List


MIN_RES_PATCH = 5
NEIGHBORS_DIST = 10
BACKBONE_ATOMS = ['C', 'CA', 'O', 'N']
MAX_STRETCH_LEN = 7
MAX_HELIX_LEN = 11


def define_motifs(pose, pdb_name, surface_only):

    if surface_only:
        surf_sel = create_layer_selector()
        selected_res = surf_sel.apply(pose)
    else:
        sel_all = rosetta.core.select.residue_selector.ResidueIndexSelector()
        for i in range(1, pose.total_residue() + 1):
            sel_all.append_index(i)
        selected_res = sel_all.apply(pose)

    include_focus = True
    neighborhood_selector = create_neighborhood_selector(NEIGHBORS_DIST, include_focus)

    motifs_created = [] # save here res indexes around which you created patches
    motifs_with_chains = [] # residues and chains [('A', '2'), ('A', '3'), ('A', '4')]

    pdbinf = pose.pdb_info()

    all_motifs_list = []
    num = 1
    for i,res in enumerate(selected_res, 1):
        motif_df = pd.DataFrame(columns=['chain_name', 'res_id','motif_name'])
        # cbetas = motif_df[motif_df.atom_full_name == 'CB'].reset_index(drop=True)
        if i <=2:
            continue
        if res:
            if i in motifs_created: # skip every second patch
                continue
            focus = create_index_selector([i])
            neighborhood_selector.set_focus(focus.apply(pose))
            atoms = pyrosetta.rosetta.utility.vector1_std_string()

            atoms.append('CA')

            # if pose.residue(i).name() == 'GLY':
            #     atoms.append('CA')
            # else:
            #     atoms.append('CB')

            neighborhood_selector.set_atom_names_for_distance_measure(atoms)
            new_motif_sel = neighborhood_selector.apply(pose)

            motif = list(map(str, core.select.residue_selector.selection_positions(new_motif_sel)))
            refine_motif(motif, pose) # it is supposed to be pointer I think
            motif_with_chains = [(pdbinf.chain(int(resn)), str(pdbinf.number(int(resn)))) for resn in motif]

            ##################################
            motifs_created.append(i+1)
            if motif_with_chains not in motifs_with_chains:
                motifs_with_chains.append(motif_with_chains)
                motif_df.chain_name = [res[0] for res in motifs_with_chains[-1]]
                motif_df.res_id = [res[1] for res in motifs_with_chains[-1]]
                motif_name = '%03d_'%num+pdb_name
                motif_df.motif_name = motif_name
                num += 1

                write_to_pdb(motif, motif_name, pose)

        all_motifs_list.append(motif_df)

    all_motifs = pd.concat(all_motifs_list).reset_index(drop=True)
    return all_motifs


def write_to_pdb(motif, motif_name, pose):
    motif_pose = Pose()
    slice_res = pyrosetta.rosetta.utility.vector1_unsigned_long()
    for res in motif:
        slice_res.append(int(res))
    core.pose.pdbslice(motif_pose, pose, slice_res)
    motif_pose.dump_pdb(motif_name + '.pdb')


def refine_motif(motif, pose):
    """Stretches are not longer than 5aa, no single residues"""
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
                res_to_remove.extend(motif[idx - (stretch_len - max_length) + 1:idx + 1]) # cut sequences longer than 5
                # stretch = stretch[:max_length]
            stretch_len = 1
            continue
        else:
            res_to_remove.append(resid)
    if stretch_len <= 1:
        res_to_remove.append(resid)
    elif stretch_len > MAX_STRETCH_LEN and stretch_len == len(motif): # cut sequences longer than 5
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


############## PyRosetta selectors from utils

def create_layer_selector():
    lay_sel = rosetta.core.select.residue_selector.LayerSelector()
    lay_sel.set_layers(0, 0, 1)
    lay_sel.set_use_sc_neighbors(0)
    lay_sel.set_ball_radius(1.35)
    # lay_sel.set_cutoffs(5.2, 2)
    return lay_sel


def create_index_selector(res_nums: List[int]):
    idx_selector = rosetta.core.select.residue_selector.ResidueIndexSelector()
    for res in res_nums:
        idx_selector.append_index(res)
    return idx_selector


def create_neighborhood_selector(cutoff, include_focus):

    neighborhood_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    neighborhood_selector.set_include_focus_in_subset(include_focus)
    neighborhood_selector.set_distance(cutoff)

    return neighborhood_selector

##############################################


def parse_pdb_file(pdb_file, motifs):
    with open(pdb_file, 'r') as f:
        lines = f.readlines()
    ids = []
    atom_full_names = []
    chain_names = []
    res_id = []
    xs = []
    ys = []
    zs = []
    for line in lines:
        if len(line.split()) > 1 and (line.startswith('ATOM') or line.startswith('HETATM')):
            ############# only backbone
            if line[12:16].strip() in BACKBONE_ATOMS:
                # if line[12:16].strip() == 'CA':
                ids.append(line[6:11].strip())
                atom_full_names.append(line[12:16].strip())
                chain_names.append(line[21])
                res_id.append(line[22:26].strip())
                xs.append(float(line[30:38].strip()))
                ys.append(float(line[38:46].strip()))
                zs.append(float(line[46:54].strip()))

    pdb_data = pd.DataFrame({'id':ids,
                             'atom_full_name':atom_full_names,
                             'chain_name':chain_names,
                             'res_id':res_id,
                             'X':xs,
                             'Y':ys,
                             'Z':zs})

    motif_data = pdb_data.merge(motifs)

    return motif_data


def split_prot(in_list, surface):
    all_prot_list = []
    if len(in_list) > 1:
        for prot in in_list:
            if os.path.isfile(prot): # This is in case that our local PDB is not updated
                try: # This one you need to research
                    motif_data = collect_motif_data(prot, surface)
                except:
                    print(prot + ' FAILED')
                    continue
                all_prot_list.append(motif_data)
            else:
                continue #########
        all_prot_data = pd.concat(all_prot_list).reset_index(drop=True)
        return all_prot_data
    prot = in_list[0]
    motif_data = collect_motif_data(prot, surface)
    return motif_data


def collect_motif_data(prot, surface):
    prot_name = os.path.splitext(os.path.basename(prot))[0]
    pose = pose_from_pdb(prot)
    motifs = define_motifs(pose,prot_name, surface)
    print("The surface was split into " + str(motifs.motif_name.nunique()) + " patches")
    motif_data = parse_pdb_file(prot, motifs)
    return motif_data


def arg_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', '-p', dest='pdb_file', default=None)
    parser.add_argument('--list', '-l', dest='pdb_list', default=None)

    return parser


def main():

    args = arg_parser().parse_args()
    pdb = args.pdb_file
    pdb_list = args.pdb_list

    if pdb_list:
        list_name = os.path.basename(pdb_list)
        with open(pdb_list, 'r') as pl:
            in_list = [line.strip() for line in pl.readlines()]
        surface = False
        outfile = list_name + '_motifs'
    elif pdb:
        in_list = [pdb]
        surface = True
        outfile = os.path.splitext(os.path.basename(pdb))[0]
    else:
        print('Provide either pdb file to split to surface patches or list of pdbs to split to motifs')
        sys.exit()

    pool = mp.Pool(mp.cpu_count())
    results = [pool.apply(split_prot, args=(in_list[i:i+50], surface)) for i in range(0, len(in_list), 50)]

    # motif_data = split_prot(in_list, surface)

    for i, motif_data in enumerate(results):
        with open(outfile + '_%s.pkl'%str(i), 'wb') as f:
            pickle.dump(motif_data, f)
    # motif_data.to_hdf(outfile, 'motifs')


if __name__ == "__main__":

    pyrosetta.init()

    main()
