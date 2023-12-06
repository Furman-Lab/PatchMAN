from pyrosetta import *
from pyrosetta.rosetta import *

import os
import sys
import pyrosetta_utils as utils
import argparse

MIN_RES_PATCH = 5 # TODO: delete, seems not to be used
NEIGHBORS_DIST = 10
BACKBONE_ATOMS = ['C', 'CA', 'O', 'N'] # TODO: delete, seems not to be used
MAX_STRETCH_LEN = 7
MAX_HELIX_LEN = 11
RATIO_ALLOWED_MASKED_RESIDUES = 0.3  # added for masking
MIN_OVERLAPPING_RESIDUES = 5  # added for masking


def create_motif_with_chains(pdbinf, motif):
	"""
	Function to create a list of tuples with the chain and the residue number
	:param pdbinf: pdb_info object for the pose
	:param motif: motif with rosetta numbering
	:return: motif as a tuple with chain and residue number
	"""
	motif_with_chains = [str(pdbinf.chain(int(resn))) + str(pdbinf.number(int(resn))) for resn in motif]
	return motif_with_chains


def check_mask_on_motif(pdbinf, motif, mask_res, debug=False):
	"""
	Function to check if the motif contains masked residues over a certain threshold
	:param motif: list of residues that are in the motif
	:param mask_res: list of residues that are masked
	:return: the motif with the chains, which contain nothing if the motif is masked
	"""
	motif_with_chains = []
	motif_length = len(motif)
	min_unmasked_length = motif_length * (1 - RATIO_ALLOWED_MASKED_RESIDUES)
	motif_pdb = utils.convert_rosetta_numbering_to_pdb_numbering(pdbinf, motif)

	# get the residues that are not in the mask
	unmasked_motif = list(set(motif_pdb) - set(mask_res))
	unmasked_length = len(unmasked_motif)

	if debug:
		print('Comparing motif and mask_res:\nmotif: ', motif_pdb, '\nmasks: ', mask_res)
		print('Unmasked motif residues: ', unmasked_motif)

	# break the loop if too many residues are masked
	if unmasked_length > min_unmasked_length:
		motif_not_masked = True
		motif_with_chains = create_motif_with_chains(pdbinf, motif)
	else:
		motif_not_masked = False

	if debug:
		ratio_of_masked = round((motif_length - unmasked_length) / motif_length * 100, 0)
		print(
			"{unmasked}/{motiflen} =  residues in the motif are unmasked ({ratio}% masked). Motif kept: {kept}".format(
				motiflen=motif_length, unmasked=unmasked_length, ratio=ratio_of_masked, kept=motif_not_masked))
		print('#####################################')

	return motif_with_chains


def modify_selected_resi(pdbinf, selected_res, pdb_number_check_resi, focus=False, debug=False):
	"""
	Function to modify the list of selected residues. If the residue is not in focus, it is removed from the list.
	If the residue is in focus, it is added to the list.
	:param pdbinf: pdb_info object for the pose
	:param selected_res: list of selected residues
	:param check_pose: list of residues to check, as a pose
	:param focus_list: list of residues in PDB numbering format
	:param focus: if True, add residues to the list, if False, remove residues from the list
	:param debug: print info
	:return: modified list of selected residues
	"""
	for i, res in enumerate(selected_res, 1):
		pdb_number = str(pdbinf.chain(int(i)))+str(pdbinf.number(int(i)))
		if res:
			if focus:  # if we want to focus on the residues, we need to set the residues that are not in focus to False
				if pdb_number not in pdb_number_check_resi:
					selected_res[i] = False
			else:  # if we want to mask the residues, we need to set the residues that are in mask to False
				if pdb_number in pdb_number_check_resi:
					selected_res[i] = False
	
	return selected_res


def define_motifs(pose, pdb_name, check_list=None, focus=True, hotspot_mode=False, edit_reslist=False, debug=False):
	"""
	Function to extract motifs from the receptor. Motifs are defined as patches of residues that are close to each other.
	Note that mask_pose and focus_pose are mutually exclusive. If you want to mask residues, you need to set mask_pose.
	:param pose: input receptor
	:param pdb_name: basename of the receptor, extracted from the filename
	:param focus: if check_list is provided, decide if it is focus or masking mode
	:param check_list: list of residues that we would like to focus on as a binding site
	:param hotspot_mode: force to use the provided residues as a center of the patch
	:param debug: print info
	:return: list of extracted motifs
	"""

	surf_sel = utils.create_layer_selector()
	selected_res = surf_sel.apply(pose)
	pdbinf = pose.pdb_info()

	# remove residues from the selected ones, but only if this is how we deal with focusing
	# otherwise, we will remove those patches that does not have residues from the focus list
	if (check_list is not None and edit_reslist):
		selected_res = modify_selected_resi(pdbinf, selected_res, check_list, focus=focus, debug=debug)

	# print out residues that are selected as focus points
	if debug:
		print_text = "Residues selected as focus points: "
		for i, res in enumerate(selected_res, 1):
			if res:
				print_text += str(pdbinf.number(i)) + " "
		print(print_text)

	include_focus_resi = True  # Alisa's original , now set to true
	neighborhood_selector = utils.create_neighborhood_selector(NEIGHBORS_DIST, include_focus_resi)
	motifs_created = []  # save here res indexes around which patches were created
	motifs_with_chains = []  # residues and chains [('A', '2'), ('A', '3'), ('A', '4')]

	all_motifs_list = []
	num = 1
	focus_res_motif = []
	for i, res in enumerate(selected_res, 1):
		# skip  2 res (too small, overlapping patches)
		# but for hotspot mode, this might delete important parts, so keep them
		if i <= 2 and not focus and not hotspot_mode:
			continue
		if res:
			# skip every second patch (to prevent too overlapping patches)
			# in hotspot mode, we force the residues in the check_list to be the center of the patch
			if i in motifs_created and ((check_list is None and focus) or not hotspot_mode):
				continue

			# For each center residue, create a patch of residues around it with the given distance
			center_res = utils.create_index_selector([i])
			neighborhood_selector.set_focus(center_res.apply(pose))
			atoms = pyrosetta.rosetta.utility.vector1_std_string()
			atoms.append('CA')
			neighborhood_selector.set_atom_names_for_distance_measure(atoms)
			new_motif_sel = neighborhood_selector.apply(pose)

			motif_residues = core.select.residue_selector.selection_positions(new_motif_sel)
			motif = list(map(str, motif_residues))
			refine_motif(motif, pose)

			motif_not_masked = True
			number_of_masked_residues = 0
			if check_list is not None and not focus:
				motif_with_chains = check_mask_on_motif(pdbinf, motif_residues, check_list, debug=debug)
			else:
				motif_with_chains = create_motif_with_chains(pdbinf, motif)
			
			# if focus mode, check if motif_with_chains have any residues in check_list
			if check_list is not None and focus:
				common_residues = set(motif_with_chains).intersection(set(check_list))
				
				# in hotspot mode, we keep every patch that overlaps with even one residue from the check_list
				# in focus mode, we keep every patch that overlaps with at least MIN_OVERLAPPING_RESIDUES residues from the check_list
				if (len(common_residues) < MIN_OVERLAPPING_RESIDUES and not hotspot_mode) or (len(common_residues) == 0 and hotspot_mode):
					if debug:
						print('WARNING: Motif does not contain any residues from the focus list, skipping')
						print('motif residues:', set(motif_with_chains))
						print('check list:', set(check_list))
					continue
			
			if len(motif_with_chains) > 0:  # add checking if motif contained masked residues
				motifs_created.append(i + 1)
				if motif_with_chains not in motifs_with_chains:
					motifs_with_chains.append(motif_with_chains)
					motif_name = '%03d_' % num + pdb_name
					num += 1

					write_to_pdb(motif, motif_name, pose)
					all_motifs_list.append(motif_name)
					focus_res_motif.append(['%03d' % num, str(pdbinf.number(int(i))), ','.join(motif)])

	with open(pdb_name + '_focus_res_motif.txt', 'w') as f:
		for l in focus_res_motif:
			f.write(pdb_name + ' ' + ' '.join(l) + '\n')

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
	# pose.display_secstruct()
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
				res_to_remove.extend(motif[idx - (stretch_len - max_length) + 1:idx + 1])  # cut long sequences
			stretch_len = 1
			continue
		else:
			res_to_remove.append(resid)
	if stretch_len <= 1:
		res_to_remove.append(resid)
	elif stretch_len > MAX_STRETCH_LEN and stretch_len == len(motif):  # cut long sequences
		if is_helix(motif[idx + 1 - stretch_len:idx + 1], s):
			max_length = MAX_HELIX_LEN
		else:
			max_length = MAX_STRETCH_LEN
		res_to_remove = motif[max_length:]
	for r in res_to_remove:
		motif.remove(r)


def is_helix(stretch, s):
	"""Check if stretch is helix"""
	if s[int(stretch[0]):int(stretch[-1]) + 1].count('H') >= 3:
		return True
	else:
		return False


def split_to_motifs():
	pyrosetta.init('-mute core.select.residue_selector.NeighborhoodResidueSelector -mute core')

	# cd into working directory - important for singularity
	os.chdir(os.getenv('work_dir'))

	# set up 3 arguments with argparse: -i, -f and -m
	parser = argparse.ArgumentParser(description = 'Create motifs from a protein structure')
	parser.add_argument('-i', '--input', help='Input PDB file', required=True)
	parser.add_argument('-f', '--focus', help='Switching between focus and masking mode', action='store_true', default=False)
	parser.add_argument('-l', '--resi_list', help='List of focus residues that should or should not form the interface, separated by comma: "A31,A56,A12"',
	                    required=False, default=None)
	parser.add_argument('-s', '--mask_focus', help='PDB file, with residues that should or should not be on the interface',
	                    required=False, default=None)
	parser.add_argument('-o', '--hotspot_mode',
	                    help='Only with focus mode. If True, all residues in the focus patch will be used as a center of patch.',
	                    required=False, action='store_true')
	parser.add_argument('-v', '--verbose', help='For debugging purposes', required=False, action='store_true',
	                    default=False)
	args = parser.parse_args()
	
	# if hospot mode is true, focus mode must be true
	if args.hotspot_mode:
		args.focus = True
	
	# load receptor pose
	pose, prot_name = utils.load_and_clean_pdb(args.input, return_name=True)
	print(args)
	
	# process inputs for mask and focus residues
	check_list = utils.parse_mask_and_focus(args)

	# split the surface into motifs
	motifs = define_motifs(pose, prot_name, check_list=check_list, focus=args.focus, hotspot_mode=args.hotspot_mode, debug=args.verbose)

	if args.hotspot_mode and len(motifs) != len(check_list) and args.focus:
		print("WARNING: Number of motifs does not match number of focus residues in hotspot mode! Do all the provided residues exist on the surface?")

	print("The surface was split into " + str(len(motifs)) + " patches")


if __name__ == "__main__":
	pyrosetta.init('-mute core.select.residue_selector.NeighborhoodResidueSelector -mute core')

	split_to_motifs()
