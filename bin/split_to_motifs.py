from pyrosetta import *
from pyrosetta.rosetta import *

import os
import sys
import argparse
import pyrosetta_utils as utils
import protocol_utils
from types import SimpleNamespace

NEIGHBORS_DIST = 10
MAX_STRETCH_LEN = 7
MAX_HELIX_LEN = 11
RATIO_ALLOWED_MASKED_RESIDUES = 0.3  # added for masking


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


def create_motif_around_center(center_res_i, namespace, debug=False):
	"""
	Function to create a motif around the center residue with the given distance
	:param pose: input pose
	:param center_res: center residue
	:param pdbinf: pdb_info object for the pose
	:param debug: print info
	:return: list of residues that are in the motif
	"""
	# Define center residue of patch
	center_res = utils.create_index_selector([center_res_i])
	namespace.neighborhood_selector.set_focus(center_res.apply(namespace.pose))

	# Use CA for selecting the residues in certain distance
	atoms = pyrosetta.rosetta.utility.vector1_std_string()
	atoms.append('CA')
	namespace.neighborhood_selector.set_atom_names_for_distance_measure(atoms)
	new_motif_sel = namespace.neighborhood_selector.apply(namespace.pose)

	# Get the residues that are in the motif
	motif_residues = core.select.residue_selector.selection_positions(new_motif_sel)
	motif = list(map(str, motif_residues))

	# Refine the motif
	refine_motif(motif, namespace.pose)

	if debug:
		print('Motif residues:', motif)

	return motif, motif_residues


def write_center_residues(namespace):
	"""Write the central residues of the motifs to a file"""
	with open(namespace.pdb_name + '_central_res_motif.txt', 'w') as f:
		for l in namespace.central_res_motif:
			f.write(namespace.pdb_name + ' ' + ' '.join(l) + '\n')


def finalize_motif(motif, motif_with_chains, namespace, i):
	"""Finalize the motif and save it to the lists"""
	namespace.motifs_created.append(i + 1)
	if motif_with_chains not in namespace.motifs_with_chains:
		namespace.motifs_with_chains.append(motif_with_chains)
		motif_name = '%03d_' % namespace.num + namespace.pdb_name
		namespace.all_motifs_list.append(motif_name)

		write_to_pdb(motif, motif_name, namespace.pose)
		namespace.central_res_motif.append(['%03d' % namespace.num, str(namespace.pdbinf.number(int(i))), ','.join(motif)])
		namespace.num += 1


def normal_mode(namespace, debug=False):
	"""
	Normal mode of the script, when we are not in masking or focusing mode
	:param debug: print info

	:return: list of extracted motifs
	"""
	for i, res in enumerate(namespace.selected_res, 1):
		if res:
			# skip every second residue to be considered as patch center (to prevent too overlapping patches)
			if i in namespace.motifs_created or i <= 2:
				continue
			# create motif around the center residue
			motif, _ = create_motif_around_center(i, namespace, debug=debug)
			motif_with_chains = create_motif_with_chains(namespace.pdbinf, motif)

			finalize_motif(motif, motif_with_chains, namespace, i)


def focus_run(check_list, namespace, hotspot=False, debug=False):
	"""
	Focus mode of the script, when we are defining a larger binding site on the surface
	:param debug: print info

	:return: list of extracted motifs
	"""
	for i, res in enumerate(namespace.selected_res, 1):
		if res:
			# skip every second residue to be considered as patch center (to prevent too overlapping patches)
			if i in namespace.motifs_created or i <= 2:
				continue

			# create motif around the center residue
			motif, _ = create_motif_around_center(i, namespace, debug=debug)
			motif_with_chains = create_motif_with_chains(namespace.pdbinf, motif)

			if len(motif_with_chains) > 0:
				common_residues = set(motif_with_chains).intersection(set(check_list))

				# in focus mode, we keep every patch that overlaps with at least 5 residues from the check_list
				if len(common_residues) > 4 or (hotspot and len(common_residues) > 2):
					finalize_motif(motif, motif_with_chains, namespace, i)
				else:
					if debug:
						print('WARNING: Motif does not contain at least 5 residues from the focus list, skipping')
						print('motif residues:', set(motif_with_chains))
						print('focus residue list:', set(check_list))



def mask_run(check_list, namespace, debug=False):
	"""
	Masking mode of the script, when we are excluding residues as a binding site on the surface
	:param debug: print info

	:return: list of extracted motifs
	"""
	for i, res in enumerate(namespace.selected_res, 1):
		if res:
			# skip every second residue to be considered as patch center (to prevent too overlapping patches)
			if i in namespace.motifs_created or i <= 2:
				continue

			# create motif around the center residue
			motif, motif_residues = create_motif_around_center(i, namespace, debug=debug)
			motif_with_chains = check_mask_on_motif(namespace.pdbinf, motif_residues, check_list, debug=debug)

			if len(motif_with_chains) > 0:
				finalize_motif(motif, motif_with_chains, namespace, i)
			else:
				if debug:
					print('WARNING: Motif contains too many residues from the mask list, skipping')
					print('motif residues:', set(motif_with_chains))
					print('mask residue list:', set(check_list))


def define_motifs(pose, pdb_name, check_list=None, focus=False, hotspot_mode=False, debug=False):
	"""
	Function to extract motifs from the receptor. Motifs are defined as patches of residues that are close to each other.
	Note that mask_pose and focus_pose are mutually exclusive. If you want to mask residues, you need to set mask_pose.
	:param pose: input receptor
	:param pdb_name: basename of the receptor, extracted from the filename
	:param focus: if check_list is provided, decide if it is focus or masking mode
	:param check_list: list of residues that we would like to focus on as a binding site
	:param hotspot_mode: force to use the provided residues also as centers of patches
	:param debug: print info

	:return: list of extracted motifs
	"""
	pose.display_secstruct() # recalculates secondary structure assignments with DSSP and stores it

	# Get surface residues
	surf_sel = utils.create_layer_selector()
	selected_res = surf_sel.apply(pose)

	# create pdb_info object
	pdbinf = pose.pdb_info()

	# create neighborhood selector for center residues
	include_focus_resi = True  # Alisa's original , now set to true
	neighborhood_selector = utils.create_neighborhood_selector(NEIGHBORS_DIST, include_focus_resi)

	helper_data = SimpleNamespace(
		pose=pose,
		pdbinf=pdbinf,
		neighborhood_selector=neighborhood_selector,
		selected_res=selected_res,
		pdb_name=pdb_name,
		num=1,
		motifs_created=[], # save here res indexes around which patches were created
		motifs_with_chains=[], # residues and chains [('A', '2'), ('A', '3'), ('A', '4')]
		all_motifs_list=[], # list of all extracted and accepted motifs
		central_res_motif=[], # central residues of the motifs
	)

	if not focus and not hotspot_mode and check_list is None:
		print('Normal mode')
		normal_mode(helper_data, debug=debug)
	elif focus and check_list and not hotspot_mode:
		print('Focus residues:', check_list)
		focus_run(check_list, helper_data, debug=debug)
	elif hotspot_mode and check_list:
		print('Hotspot mode:', check_list)
		focus_run(check_list, helper_data, hotspot=True, debug=debug)
	elif not focus and check_list and not hotspot_mode:
		print('Mask residues:', check_list)
		mask_run(check_list, helper_data, debug=debug)
	else:
		print('ERROR: Wrong combination of arguments: focus, hotspot_mode, check_list')
		sys.exit(1)

	write_center_residues(helper_data)

	return helper_data.all_motifs_list


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

	# load receptor pose
	pose, prot_name = utils.load_and_clean_pdb(args.input, return_name=True)
	print(args)

	# process inputs for mask and focus residues
	check_list = utils.parse_mask_and_focus(args)

	# split the surface into motifs
	motifs = define_motifs(pose, prot_name, check_list=check_list, focus=args.focus, hotspot_mode=args.hotspot_mode, debug=args.verbose)
	print("The surface was split into " + str(len(set(motifs))) + " patches")

	# Listing motif files and creating motif_list
	protocol_utils.create_motif_list(os.path.basename(args.input)[:-4], 'motif_list')

	# Create search files for MASTER
	protocol_utils.run_createPDS('motif_list')


if __name__ == "__main__":
	pyrosetta.init('-mute core.select.residue_selector.NeighborhoodResidueSelector -mute core')

	split_to_motifs()
