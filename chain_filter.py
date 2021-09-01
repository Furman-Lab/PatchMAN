#!/vol/ek/Home/alisa/python3.5/bin/python3.5

from sklearn.neighbors import NearestNeighbors
import numpy as np
from pandas import DataFrame

BACKBONE_ATOMS = ['N','O','C','CA']

CLOSEST_DISTANCE_PERCENTILE = 0.45
MAX_NON_INT_RESIDUES = 3


def parse_pdb_file(patch_file):
    with open(patch_file, 'r') as f:
        lines = f.readlines()
    atom_full_names = []
    chain_names = []
    xs = []
    ys = []
    zs = []
    res_ids = []
    for line in lines:
        if line.startswith('ATOM') or line.startswith('HETATM'):
            atom_full_names.append(line[12:16].strip())
            chain_names.append(line[21].strip())
            res_ids.append(int(line[22:26].strip()))
            xs.append(float(line[30:38].strip()))
            ys.append(float(line[38:46].strip()))
            zs.append(float(line[46:54].strip()))

    pdb_data = DataFrame({'atom_full_name': atom_full_names,
                          'chain_name': chain_names,
                          'res_id': res_ids,
                          'X': xs,
                          'Y': ys,
                          'Z': zs})
    return pdb_data


def find_relevant_chains(prot_complex_file, clash_dist, interaction_dist):
    prot_data = parse_pdb_file(prot_complex_file)
    receptor = prot_data[prot_data.chain_name == 'A'].copy()

    pep_template = prot_data[prot_data.chain_name == 'B'].copy()
    # assert pep_template.chain_name.nunique() == 1, '{}, {}: More than 1 chain_name in template" {}'.format(
    #     prot_complex_file, receptor_chain, pep_template.chain_name.unique()) # shouldn't happen. checking earlier

    all_names = prot_data.atom_full_name.unique()
    h_names = [name for name in all_names if 'H' in name and name != 'OH'] # all the hydrogens

    protein_backbone = receptor[receptor.atom_full_name.isin(BACKBONE_ATOMS)]
    protein_backbone_coords = protein_backbone[['X','Y','Z']].values

    receptor_nonh_atoms = receptor[~receptor.atom_full_name.isin(h_names)].copy()
    receptor_nonh_atom_coords = receptor_nonh_atoms[['X','Y','Z']].values

    template_backbone = pep_template[pep_template.atom_full_name.isin(BACKBONE_ATOMS)].copy()
    template_backbone_coords = template_backbone[['X','Y','Z']].values

    if template_backbone.empty:
        print('No peptide template for threading')
        return False

    template_nonh_atoms = pep_template[~pep_template.atom_full_name.isin(h_names)].copy()
    template_nonh_atom_coords = template_nonh_atoms[['X','Y','Z']].values

    print('total atoms: protein = {}, template = {}'.format(len(receptor), len(pep_template)))
    print('backbone atoms: protein = {}, template = {}'.format(len(protein_backbone), len(template_backbone)))

    sc_neighbors = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(receptor_nonh_atom_coords)
    closest_protein_atoms = sc_neighbors.kneighbors(template_nonh_atom_coords) # for interaction screening

    bb_neighbors = NearestNeighbors(n_neighbors=1, algorithm='kd_tree').fit(protein_backbone_coords)
    closest_bb_atoms = bb_neighbors.kneighbors(template_backbone_coords) # for clashes

    template_nonh_atoms['closest_dist'] = closest_protein_atoms[0][:, 0]
    template_backbone['closest_dist'] = closest_bb_atoms[0][:, 0]

    # Check for gaps in the peptide. Assume that CA-CA distance in a peptide bond is < 4A
    peptide_CAs = template_backbone[template_backbone.atom_full_name == 'CA']
    for i, coord in enumerate(peptide_CAs[['X','Y','Z']].values):
        if i < len(peptide_CAs[['X','Y','Z']].values) - 1:
            if np.linalg.norm(peptide_CAs[['X','Y','Z']].values[i+1] - peptide_CAs[['X','Y','Z']].values[i]) > 4:
                print("The peptide is broken")
                return False

    # Check for clash
    min_bb_dist = np.round(template_backbone.closest_dist.min(), 4)
    if min_bb_dist < clash_dist:
        print('CLASH detected! min distance is {} (< clash_distance={})'.format(min_bb_dist, clash_dist))
        return False
    else:
        # Check for interaction
        num_template_residues = template_nonh_atoms.res_id.nunique()
        res_min_dists = template_nonh_atoms.groupby('res_id').closest_dist.min().reset_index()

        res_min_dists = res_min_dists.sort_values('res_id')

        num_non_interacting_residues = len(res_min_dists[res_min_dists.closest_dist > interaction_dist])
        non_interacting_pct = num_non_interacting_residues / float(num_template_residues)
        print('{} out of {}({}%) template backbone residues have min atom distance > interaction dist={}'.format(
            num_non_interacting_residues, num_template_residues, np.round(non_interacting_pct,2), interaction_dist))

        if np.all(res_min_dists.head(2).closest_dist > interaction_dist):
            return False
        if np.all(res_min_dists.tail(2).closest_dist > interaction_dist):
            return False # two non-interacting res at the end
        return non_interacting_pct < CLOSEST_DISTANCE_PERCENTILE



