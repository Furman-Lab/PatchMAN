from pyrosetta import *
from pyrosetta.rosetta import *
import numpy as np

# This function currently returns most of the protein as surface, it is not clear why.
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


def extract_check_list_from_file(filepath):
    """
    Extract residue numbers and chain IDs from PDB or CIF files.

    Args:
        filepath (str): Path to the PDB or CIF file

    Returns:
        str: Comma-separated list of chain ID + residue number (e.g., 'A13,A34,A56')

    Raises:
        ValueError: If required column information cannot be found in the CIF file header
    """
    # Determine file type from extension
    file_type = filepath.lower().split('.')[-1]
    
    residues = set()  # Use a set to avoid duplicates
    
    if file_type == 'pdb':
        with open(filepath, 'r') as file:
            for line in file:
                if line.startswith('ATOM'):
                    # PDB format: columns are fixed width
                    # Chain ID is typically at position 21
                    # Residue number is typically at positions 22-26
                    chain_id = line[21].strip()
                    res_num = line[22:26].strip()
                    residues.add(f"{chain_id}{res_num}")
    
    elif file_type == 'cif':
        # For CIF files, we need to find the column indices from the header
        auth_asym_id_col = -1  # Chain ID column
        auth_seq_id_col = -1  # Residue number column
        
        with open(filepath, 'r') as file:
            header_lines = []
            in_header = False
            in_atom_data = False
            
            for line in file:
                line = line.strip()
                
                # Start collecting header when we see atom_site fields
                if line.startswith('_atom_site.'):
                    in_header = True
                    header_lines.append(line.split('.')[-1])
                
                # When we exit the header section, process the collected headers
                elif in_header and not line.startswith('_'):
                    # Process header to find column indices
                    chain_col = header_lines.index('auth_asym_id') if 'auth_asym_id' in header_lines else -1
                    res_num_col = header_lines.index('auth_seq_id') if 'auth_seq_id' in header_lines else -1
                    
                    # If they are still -1, search for another key
                    if chain_col == -1:
                        chain_col = header_lines.index('label_asym_id') if 'label_asym_id' in header_lines else -1
                    if res_num_col == -1:
                        res_num_col = header_lines.index('label_seq_id') if 'label_seq_id' in header_lines else -1
                    
                    # Raise error if we couldn't find the required columns
                    if chain_col == -1 or res_num_col == -1:
                        raise ValueError(
                            "Could not find required column information in CIF header. "
                            "Need 'auth_asym_id'/'label_asym_id' for chain ID and "
                            "'auth_seq_id'/'label_seq_id' for residue number."
                        )
                    
                    in_header = False
                    in_atom_data = True
                
                # Process ATOM records
                if in_atom_data and line.startswith('ATOM') and 'CA' in line:
                    columns = line.split()
                    
                    # Make sure we have enough columns
                    if len(columns) > max(chain_col, res_num_col):
                        chain_id = columns[chain_col]
                        res_num = columns[res_num_col]
                        residues.add(f"{chain_id}{res_num}")
                    else:
                        # This is an unexpected format error in the ATOM record
                        raise ValueError(f"Warning: Skipping malformed ATOM record: {line}")
    
    return residues


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
        check_list = extract_check_list_from_file(args.mask_focus)

    print('Mask/focus residues are: ' + ','.join(check_list))

    return check_list

def load_and_clean_pdb(pdb_file, return_name=False):
    """Load and clean pdb file
    :param pdb_file: pdb file
    :param return_name: return name of
    :return: pose and name of  if asked for
    """
    # TODO: mask_focus crashes on first loading for some reason, even if a valid cif file is provided

    # print working directory
    print("cleaning: ", pdb_file)
    prot_name = os.path.splitext(os.path.basename(pdb_file))[0]

    # this is because cleanATOM doesnt correctly handle cif files
    if pdb_file.endswith('.cif'):
        pose = pyrosetta.rosetta.core.import_pose.pose_from_file(pdb_file)
        pdb_file = f'{prot_name}.pdb'
        pose.dump_pdb(pdb_file)

    toolbox.cleaning.cleanATOM(pdb_file)
    clean_filename = f'{prot_name}.clean.pdb'

    print(f'Loading : {clean_filename}')
    pose = core.import_pose.pose_from_file(clean_filename)

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
