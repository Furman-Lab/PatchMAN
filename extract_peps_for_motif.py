#!/vol/ek/Home/alisa/python3.5/bin/python3

from pyrosetta import *
from pyrosetta.rosetta import *
import sys
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.SubsMat import MatrixInfo as matlist
import time
from os import system
from os import path
import argparse

sys.path.insert(1, '/vol/ek/Home/alisa/scripts/motifs/motif_scan')

import pyrosetta_utils as utils
from chain_filter import find_relevant_chains
from rosetta_protocols import fixbb_design

##### DEBUG
GAP_ALLOWED = 3
ELONGATE_BY_MAX = 4
NEIGHBORS_DIST = 8
INTERACTION_DIST = 5
CLASH_DIST = 2

OVERALL_MATCHES = 100


"""Receives a list of matches for *1* motif, the receptor pdb file and the peptide sequence.
Create initial complexes: extract peptides from proteins with motifs similar to the query, thread the pepseq with fixbb.

*Take first 50 matches --> low RMSD, and last 50 matches --> depends on the RMSD cutoff"""


def extract_templates_for_motif(matches, pepseq, plen, patch, receptor_pose, scrfxn, design):
    """For each motif there are N matches. Currently I limit them to the 1000 best RMSD matches. Probably should
    sample more distant matches too."""
    single_motif_complexes = 0

    start_motif = time.time()
    print('Begin generating complexes for a motif')

    patch_pose = pose_from_pdb(patch)

    patch_indices = [patch_pose.pdb_info().pose2pdb(i).split()[0] for i in range(1, patch_pose.total_residue() + 1)]

    patch_name = patch.split('_')[0]
    log_name = patch_name + '.log'

    with open(log_name, 'w') as log:
        log.write('Patch: %s\n'%patch_name)

    print(matches)
    for match in matches:
        rmsd = match[0]
        motif_pdb = match[1][-4:]
        motif_stretches = match[2]
        t = match[3]
        R = match[4]

        indices = [r + 1 for m_stretch in motif_stretches for r in m_stretch]  # the numbering in master output is from 0

        match_to_report = motif_pdb + ': ' + ','.join([str(i) for i in indices])

        if not path.isfile(motif_pdb.upper() + '.clean.pdb'):
            pdb_pose = toolbox.pose_from_rcsb(motif_pdb)
        else:
            pdb_pose = pose_from_pdb(motif_pdb.upper() + '.clean.pdb')

        chain_breaks = []
        for jump in range(1, pdb_pose.num_chains()):
            chain_breaks.append(pdb_pose.chain_end(jump))

        try:
            env_res, env_res_with_chain, tot_res, pdb_pose = res_around_patch(indices, pdb_pose)
        except RuntimeError:  # why?
            continue

        # filter short stretches
        stretches = choose_stretches_only(env_res, chain_breaks, plen, tot_res)  # elongate and filter afterwards

        if not stretches:
            continue

        superimposed_pose = superimpose_using_RT(t, R, pdb_pose)
        for i, stretch in enumerate(stretches):
            complex_name = patch_name + '_' + motif_pdb + '_' + str(indices[0]) + '_%02d' % i + '.pdb'
            complex_pose = create_complex(receptor_pose, superimposed_pose, stretch, complex_name, patch_indices)

            if complex_pose:

                if not design:
                    if thread_pepseq(complex_name, complex_pose, pepseq, scrfxn):
                        single_motif_complexes += 1

                    system('rm -f {}'.format(complex_name))

                pep_template_seq = complex_pose.chain_sequence(2)
                motif_seq, patch_seq = compare_motif_seq_id(patch_pose, pdb_pose, indices)
                complex_inf = print_inf(complex_name, match_to_report, motif_seq, patch_name, patch_seq,
                                        pep_template_seq, pepseq, stretch,rmsd)  # print the information about the motif and patch + alignments of the motifs and peps
                with open(log_name, 'a') as log:
                    log.write(complex_inf)
            else:
                system('rm -f {}'.format(complex_name))

    print('motif ran for %s min, %s complexes were created'%(str((time.time()-start_motif)/60),
                                                                str(single_motif_complexes)))

    return single_motif_complexes


def thread_pep(complex_name, complex_pose, indices, log_name, match_to_report, patch_name, patch_pose, pdb_pose, pepseq,
               rmsd, scrfxn, stretch):
    thread_pepseq(complex_name, complex_pose, pepseq, scrfxn)
    motif_seq, patch_seq = compare_motif_seq_id(patch_pose, pdb_pose, indices)
    pep_template_seq = complex_pose.chain_sequence(2)
    complex_inf = print_inf(complex_name, match_to_report, motif_seq, patch_name, patch_seq,
                            pep_template_seq, pepseq, stretch,
                            rmsd)  # print the information about the motif and patch + alignments of the motifs and peps
    with open(log_name, 'a') as log:
        log.write(complex_inf)


def print_inf(complex_name, match_to_report, motif_seq, patch_name, patch_seq, pep_template_seq, pepseq, stretch, rmsd):

    pep_alignments = pairwise2.align.globaldx(pepseq, pep_template_seq, matlist.blosum62)
    motif_alignments = pairwise2.align.globaldx(patch_seq, motif_seq, matlist.blosum62)

    inf = 'Complex_name = {complex}; patch = {patch}; match = {match}; stretch = {stretch}\n'.format(
        complex=complex_name, patch=patch_name, match=match_to_report, stretch=','.join([str(r) for r in stretch])) + \
          'Pep seq: %s\n' % pepseq + 'Template pep seq: %s\n' % pep_template_seq + format_alignment(*pep_alignments[0]) + \
          'Patch seq: %s\n' % patch_seq + 'Motif seq: %s\n' % motif_seq + 'Motif RMSD: %s\n' % rmsd + \
          format_alignment(*motif_alignments[0]) + \
          '====================================================\n'

    return inf


def compare_motif_seq_id(patch_pose, motif_pose, indices):
    patch_seq = patch_pose.sequence()
    motif_pdb_seq = motif_pose.sequence()
    motif_seq = ''
    for i in indices:
        try:
            motif_seq += motif_pdb_seq[i]
        except IndexError:
            print('IndexError!\n')
            print('Patch: %s\n'%patch_pose.pdb_info().name())
            print('Motif: %s\n'%motif_pose.pdb_info().name())
            return None, None
    return motif_seq, patch_seq


def thread_pepseq(complex_pose_name, complex_pose, pepseq, scrfxn):
    pep_num = list(range(complex_pose.chain_begin(2), complex_pose.chain_end(2) + 1))
    idx_sel = utils.create_index_selector(pep_num)
    pep_sele = idx_sel.apply(complex_pose)

    if not fixbb_design(pep_sele, complex_pose_name, pepseq, scrfxn) != 1:
        print('fixbb failed on complex {}'.format(complex_pose_name))


def create_complex(receptor, pose_to_cut, pep, complex_name, patch_indices):
    complex_pose = Pose()
    complex_pose.assign(receptor)

    core.pose.append_subpose_to_pose(complex_pose, pose_to_cut, int(pep[0]), int(pep[-1]), True)

    complex_pose.dump_pdb(complex_name)

    #### I am not asking for pepchain, because for some reason it returns '^'
    #### In the chain filtering step I will take the chain which is not a receptor chain as a peptide chain

    # New filtering: clashes and non-interacting chains
    is_possible_chain = find_relevant_chains(complex_name, CLASH_DIST, INTERACTION_DIST, patch_indices)

    if is_possible_chain:
        return complex_pose
    else:
        return False


def superimpose_using_RT(t, R, pdb_pose):
    """Superimpose: Multiply the matching motif first by its RT and then by the inverse RT of the surface patch"""
    new_pose = Pose()
    new_pose.assign(pdb_pose)

    R_rosetta = create_r_matrix(R)
    t_rosetta = numeric.xyzVector_double_t(t[0], t[1], t[2])

    new_pose.apply_transform_Rx_plus_v(R_rosetta, t_rosetta)
    # new_pose.dump_pdb('tmp7.pdb')

    return new_pose


def create_r_matrix(R_l):
    """create matrix for transformation (from TM-align)"""
    x1, y1, z1 = R_l[0]
    x2, y2, z2 = R_l[1]
    x3, y3, z3 = R_l[2]
    u1 = numeric.xyzVector_double_t(x1, y1, z1)
    u2 = numeric.xyzVector_double_t(x2, y2, z2)
    u3 = numeric.xyzVector_double_t(x3, y3, z3)
    R = numeric.xyzMatrix_double_t()
    R.clear()
    R.row_x(u1)
    R.row_y(u2)
    R.row_z(u3)
    return R


def res_around_patch(indices, pdb_pose):
    """Found the residues around the patch to use as a template for threading"""

    pose_pdb_inf = pdb_pose.pdb_info()

    motif_sel = utils.create_index_selector([int(idx) for idx in indices])
    include_focus = False
    env_sel = select_environment(pdb_pose, motif_sel, include_focus, NEIGHBORS_DIST)
    env_res = list(map(str, core.select.residue_selector.selection_positions(env_sel)))
    tot_res = pdb_pose.total_residue()
    env_res_with_chain = [pose_pdb_inf.pose2pdb(int(res)) for res in env_res]
    return env_res, env_res_with_chain, tot_res, pdb_pose


def select_environment(pose, selector, include_focus, cutoff):

    pose_res = selector.apply(pose)
    neighborhood_selector = utils.create_neighborhood_selector(cutoff, include_focus)

    neighborhood_selector.set_focus(pose_res)

    env = neighborhood_selector.apply(pose)

    return env


def choose_stretches_only(env_res, chain_breaks, peplen, tot_res):
    stretches = []
    stretch = [env_res[0]]
    for i, res in enumerate(env_res):
        if i == len(env_res) - 1:  # The END
            if len(stretch) < 2:
                break
            else:
                for st in elongate_stretch(stretch, peplen, tot_res, chain_breaks):
                    stretches.append(st)
        elif len(stretch) >= peplen: # the stretch is long enough
            stretches.append(stretch)
            stretch = [env_res[i + 1]]
        # check if the res are sequential and that they are belong to the same chain
        elif int(res) + 1 == int(env_res[i + 1]) \
                and int(res) + 1 not in chain_breaks:
            stretch.append(env_res[i + 1])
        elif len(stretch) >= 2:
            for st in elongate_stretch(stretch, peplen, tot_res, chain_breaks):
                stretches.append(st)
            stretch = [env_res[i + 1]]
        else:
            stretch = [env_res[i + 1]]
    return stretches


def elongate_stretch(stretch, peplen, tot_res, chain_breaks):
    """Elongate stretches that are shorter than the peptide seq"""
    dif = peplen - len(stretch)
    new_stretch = []
    stretches = []

    first_res = int(stretch[0]) - dif  # from which we CAN start??
    while first_res <= 0:
        first_res += 1
        if first_res == int(stretch[0]):
            break
    for i in range(first_res, int(stretch[0]) + 1):
        if i+peplen > tot_res + 1:
            break
        for j in range(i, i+peplen):
            new_stretch.append(j)
            if j in chain_breaks:
                break
        if len(new_stretch) == peplen:
            stretches.append(new_stretch)
        new_stretch = []

    return stretches


def parse_matches(match_file):
    """Parse MASTER output, extract patch/motif res numbers and RT matrices"""
    with open(match_file.strip(), 'r') as m:
        all_all_matches = m.readlines()

    print("%s matches were found for this motif"%str(len(all_all_matches)))

    if len(all_all_matches) > OVERALL_MATCHES:
        all_matches = all_all_matches[:int(OVERALL_MATCHES/2)] + all_all_matches[-int(OVERALL_MATCHES/2):]
    else:
        all_matches = all_all_matches

    matches = []
    for line in all_matches:
        rmsd = float(line.split()[0])
        if rmsd < 0.01: # skip identical matches
            continue
        pdb = line.split()[1].split('/')[-1][:-4]
        stretches = []
        match_ranges = line[line.find('[') + 1:line.find(']')]
        match_res = match_ranges.split(',')
        for k in range(0, len(match_res), 2):
            stretches.append(list(range(int(match_res[k].lstrip().lstrip('(')), int(match_res[k + 1].rstrip(')')))))
        t = [float(val) for val in line[line.find('T') + 2:line.find('U')].strip().split()]
        u_tmp = [float(val) for val in line[line.find('U') + 2:line.find('===')].strip().split()]
        u1, u2, u3 = [], [], []
        for j, v in enumerate(u_tmp):
            if j < 3:
                u1.append(v)
            elif j < 6:
                u2.append(v)
            else:
                u3.append(v)
        u = [u1, u2, u3]

        matches.append((rmsd, pdb, stretches, t, u))

    return matches


def arg_parser():
    parser = argparse.ArgumentParser(description='Extract templates from pdb and thread the sequence/design a pep')
    parser.add_argument('--match_list', '-m', dest='matchl', default=None)  # list of matches from master
    parser.add_argument('--design', '-d', dest='design', action='store_true', default=False)  # bool
    parser.add_argument('--receptor', '-r', dest='rec', default=None)  # receptor pdb
    parser.add_argument('--patch', '-a', dest='patch', default=None)  # patch pdb
    parser.add_argument('--peptide', '-p', dest='pep', default=None)  # peptide seq for docking
    parser.add_argument('--peplen', '-l', dest='plen', default=None)  # peptide length for design

    return parser


def main():
    args = arg_parser().parse_args()
    all_matches = parse_matches(args.matchl) # here pass the list of all matches for 1 motif
    pep = args.pep
    receptor = args.rec
    patch = args.patch


    start_all = time.time()

    scrfxn = create_score_function('ref2015')
    receptor_pose = pose_from_pdb(receptor)

    if not args.design:
        with open(pep, 'r') as peptide:
            peptide_seq = peptide.readlines()
            if peptide_seq[0][0] == '>':
                pepseq = peptide_seq[1].strip()
            else:
                pepseq = peptide_seq[0].strip()
        plen = len(pepseq)
    else:
        pepseq=None
        plen = int(args.plen) # peptide length (for design)

    extract_templates_for_motif(all_matches, pepseq, plen, patch, receptor_pose, scrfxn, args.design)

    print('All templates were generated in %s min'%(str((time.time()-start_all)/60)))


if __name__ == "__main__":

    init()  # for pyrosetta

    main()
