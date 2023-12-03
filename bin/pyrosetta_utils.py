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


def parse_mask_and_focus(args):
    """Parse mask and focus residues
    :param args: command line arguments
    :return: list of residues to mask or to focus on
    """
    
    # check for case when no residues are provided
    if args.resi_list is None and args.mask_focus is None:
        return None
    # check for invalid input
    elif args.resi_list is not None and args.mask_focus is not None:
        print('ERROR: Please provide either a list or a pdb file not both!')
        sys.exit(1)
    # process residue list
    elif args.resi_list is not None:
        check_list = args.resi_list.split(',')
    # process pdb file
    elif args.mask_focus is not None:
        pose = load_and_clean_pdb(args.mask_focus)
        check_list = convert_rosetta_numbering_to_pdb_numbering(pose.pdb_info())
    
    if not args.focus and args.hotspot_mode:
        print('WARNING: Hotspot mode does not work with masking, omitting hotspot mode')
        args.hotspot_mode = False
    
    print('mask/focus residues are: ' + ','.join(check_list))
    
    return check_list

def load_and_clean_pdb(pdb_file, return_name=False):
	"""Load and clean pdb file
	:param pdb_file: pdb file
	:param return_name: return name of
	:return: pose and name of  if asked for
	"""
	# print working directory
	toolbox.cleaning.cleanATOM(pdb_file)
	prot_name = os.path.splitext(os.path.basename(pdb_file))[0]
	print('Loading ' + prot_name + '.clean.pdb')
	pose = pose_from_pdb(prot_name + '.clean.pdb')

	if return_name:
		return pose, prot_name
	else:
		return pose


def convert_rosetta_numbering_to_pdb_numbering(pdbinf, resn_list=None):
    """
	Function to convert the numbering of the residues from Rosetta numbering to PDB numbering
	:param pose: input pose
	:param resn_list: list of residue numbers in Rosetta numbering
	:return: list of residue numbers in PDB numbering
	"""
    if not resn_list:
        resn_list = [x for x in range(1, pdbinf.nres() + 1)]
    pdb_number_resn_list = []
    
    for resn in resn_list:
        pdb_number_resn_list.append(str(pdbinf.chain(int(resn))) + str(pdbinf.number(int(resn))))
    return pdb_number_resn_list