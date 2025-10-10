import argparse
import difflib
import math
import os
from shutil import copyfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pymol
from pymol import cmd
import glob

import pyrosetta
import pyrosetta.distributed.packed_pose as packed_pose
import pyrosetta.io as io
import pyrosetta.rosetta.core as core
from pyrosetta.io.silent_file_map import SilentFileMap
from pyrosetta.rosetta.core.id import AtomID
from pyrosetta.rosetta.core.pose import Pose, tag_from_pose, tag_into_pose, get_chains
from pyrosetta.rosetta.core.scoring import superimpose_pose, gdtha
from pyrosetta.rosetta.protocols.cluster import ClusterPhilStyle
from pyrosetta.rosetta.std import map_core_id_AtomID_core_id_AtomID

import warnings
warnings.filterwarnings("ignore", category=UserWarning)


def create_fpd_scoring_mover(native_pose):
        from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser

        # print(f"Creating FlexPepDock mover")

        # FlexPepDock XML definition
        xml = f"""
                <ROSETTASCRIPTS>
                        <SCOREFXNS>
                                <ScoreFunction name="sfxn" weights="ref2015"/>
                        </SCOREFXNS>
                        <MOVERS>
                                <FlexPepDock name="fpd" pep_refine="0" extra_scoring="1" lowres_preoptimize="0" min_receptor_bb='0'/>
                        </MOVERS>
                        <PROTOCOLS>
                                <Add mover_name="fpd" />
                        </PROTOCOLS>
                </ROSETTASCRIPTS>
        """

        # Run RosettaScripts
        parser = RosettaScriptsParser()
        default_options = pyrosetta.rosetta.basic.options.process()
        tag = parser.create_tag_from_xml_string(xml, default_options)
        fpd_protocol = parser.parse_protocol_tag(tag, default_options).get_mover(1)
        fpd_protocol.set_default()  # this will take care of unboundrot

        fpd_protocol.set_native_pose(native_pose)

        return fpd_protocol

def rescore_models(silent_file, scores_sorted, tags=None, output_file='rescore.sc'):
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job
    # Remove old RMSD columns
    rms_columns = [col for col in scores_sorted.columns if col.startswith('rms')]
    scores_filtered = scores_sorted.drop(columns=rms_columns)

    # Prepare top scoring decoy
    top_scoring_tag = scores_sorted.head(1).description.iloc[0]
    top_scoring_pose = read_poses_from_silentfile(silent_file, tags=[top_scoring_tag], sort_by_score='reweighted_sc')[0]

    poses = read_poses_from_silentfile(silent_file)
    fpd_scoring_mover = create_fpd_scoring_mover(top_scoring_pose)

    #Rescore
    scores_to_save = ['rmsALL', 'rmsALL_allIF', 'rmsALL_if', 'rmsBB', 'rmsBB_allIF', 'rmsBB_if', 
                      'rmsSC_allIF', 'rmsPHIPSI', 'rmsPHIPSI', 'rmsPHIPSI_if', 'rmsCA', 'rmsCA_if']
    more_scores_to_save = ['score', 'fa_atr', 'fa_rep', 'fa_sol', 'fa_dun', 'hbond_sc']
    scores = []
    all_poses = []
    for pose in poses:
        tag = pyrosetta.rosetta.core.pose.tag_from_pose(pose)
        
        # before = pose.scores['rmsBB_if'] # debug
        fpd_scoring_mover.apply(pose)
        extra_scores = dict(get_string_real_pairs_from_current_job())
        extra_scores['description'] = tag
                
        # we also need them in the poses we work with, otherwise the final results will contain the old values
        if tag in tags:
            for k, v in extra_scores.items():
                if not k.startswith('best'):
                    pose.scores[k] = v
            all_poses.append(pose)

        # we want to also provide these, but these are not from FPD
        for score_type in more_scores_to_save:
            extra_scores[score_type] = pose.scores[score_type]
        scores.append(extra_scores)

    # Save calculated RMSD-s
    rescore_df = pd.DataFrame(scores).sort_values(by='reweighted_sc', ascending=True)
    rescore_df.to_csv(output_file, sep='\t', index=False, header=False)

    return rescore_df, all_poses

def filter_and_select_top_structures(scores, tags_to_remove=[], score_type='reweighted_sc', topX=None):
    """
    Filters the score file to exclude certain descriptions and selects top X percent structures based on score.

    :param score_file: Path to the score file (e.g., "sorted.sc").
    :param percent: Percentage of top structures to select (default is 1%).
    :return: A DataFrame containing the filtered top X percent structures.
    """
    # Filter possible disulfide bridges or where description matches additional removal criteria
    df_filtered = scores[~(
            (scores['reweighted_sc'] < -1000000) &
            (~scores['description'].isin(tags_to_remove))
    )]

    # Calculate the number of rows to take for the top X percent
    num_lines = len(df_filtered)
    if not topX:
        topX = 0.01

    if 0 < topX <= 1:  # if it is in precentage
        num_str = int(num_lines * topX)
    elif topX <= 0:
        raise ValueError('Invalid value for selecting the top X for clustering: needs to be > 0')
    else:
        num_str = topX
    num_str = max(num_str, 50) # but always cluster at least 50
    print(f'Selecting {num_str} structures for clustering')
    # if num_str > num_lines, it will just return all in pandas

    # Select the top X based on sorted scores
    top_percent_df = df_filtered.nsmallest(num_str, score_type)
    top_percent_df.to_csv('top1percent', sep='\t')

    return top_percent_df


def read_pose(silent_struct, tag, score=None, sort_by_score='reweighted_sc'):
    pose = Pose()
    silent_struct.fill_pose(pose)
    pyrosetta.rosetta.core.pose.tag_into_pose(pose, tag)
    silent_struct.energies_into_pose(pose)

    # This is for making our life easier when sorting later. It will not be written out.
    silent_struct.energies_into_pose(pose)
    if score: # for non-rosetta inputs
        pose.scores['silent_score'] = score
    elif sort_by_score in pose.scores.keys():
        pose.scores['silent_score'] = pose.scores[sort_by_score]
    else:
        print (f'{sort_by_score} cannot be read from silent file, defaulting to total_score')

    return pose


def read_poses_from_silentfile(silent_file, tags=None, sort_by_score='reweighted_sc'):
    """
    :param silent_file: the silent file to read
    :param tags: list of tags to extract. If not defined, will use all the decoys
    :param sort_by_score:

    :return:
    """
    all_poses = []

    # Setup reading from silent file
    silent_file_map = SilentFileMap(silent_file)

    if tags:
        for tag in tags:
            # Read the structure from the file by tag
            silent_struct = silent_file_map.get(tag)
            pose = read_pose(silent_struct, tag)
            all_poses.append(pose)
    else:
        for tag, silent_struct in silent_file_map.items():
            pose = read_pose(silent_struct, tag)
            all_poses.append(pose)

    return all_poses


def read_poses_from_list(list_of_pdbs, list_of_scores=None, start_pep=None, end_pep=None, add_rec=True, invert_score=False):
    all_poses = []

    if type(list_of_pdbs) == dict:
        dict_of_inputs = list_of_pdbs
    elif type(list_of_pdbs) == list and type(list_of_scores) == list:
        if len(list_of_pdbs) != len(list_of_scores):
            raise ValueError('List of PDBs and scores must be of the same length')
        dict_of_inputs = dict(zip(list_of_pdbs, list_of_scores))


    for i, pose_path in enumerate(list_of_pdbs):
        pose = io.pose_from_file(pose_path)
        
        # if we need to cut the peptide, put it in a new pose
        if start_pep and end_pep:
            start_pep_pose = pose.pdb_info().pdb2pose('B', start_pep)
            end_pep_pose = pose.pdb_info().pdb2pose('B', end_pep)
            if end_pep_pose == 0:
                raise ValueError(f'A peptide residue with number {end_pep} does not exist in {pose_path}.')
            # build the new pose
            if add_rec:
                pose_new = pose.split_by_chain(1)
                core.pose.append_subpose_to_pose(pose_new, pose, start_pep_pose, end_pep_pose, True)
            else:
                pose_new = Pose()
                core.pose.append_subpose_to_pose(pose_new, pose, start_pep_pose, end_pep_pose)
            pose = pose_new
            del pose_new # just to be on the safe side
        elif not add_rec:
            pose_new = pose.split_by_chain(2) # add only the full peptide
            pose = pose_new
            del pose_new # just to be on the safe side
        else:
            pass # this will use the whole complex for clustering
            
        
        # we need to name the pose
        tag = os.path.basename(pose_path).split('.')[0]
        tag_into_pose(pose, tag)
        
        # if a score was input, store it in the pose
        if list_of_scores:
            score = list_of_scores[i]
            # inverting is for cases when the largest value is better
            pose.scores['silent_score'] = -score if invert_score else score

        all_poses.append(pose)
    
    return all_poses

def extract_peptides(list_of_poses):
    peptides = []
    for pose in list_of_poses:
        peptides.append(pose.split_by_chain(2))
    return peptides


def extract_and_sort_score_file(score_file='score.sc', output_file="sorted.sc"):
    """
    Extract columns from the score file (score.sc), sort by 'reweighted_sc', and return the result as a pandas DataFrame.

    :param score_file: Path to the score file (e.g., "score.sc").
    :param output_file: Path to the output file (e.g., "sorted.sc"). This file will be created only if not already exists.
    :return: pandas DataFrame containing the sorted data.
    """
    # List of scores to extract
    scores_types = ['reweighted_sc', 'I_sc', 'rmsBB', 'rmsBB_if', 'description']
    
    # check if score_file starts with SEQUENCE, leave out the first line
    f = open(score_file)
    first = f.readline()
    skiprows = 1 if first.startswith('SEQUENCE:') else 0
    f.close()

    scores = pd.read_csv(score_file, sep='\s+', header=0, skiprows=skiprows, engine='python')
    scores = scores[~scores['SCORE:'].str.contains('SEQUENCE:', case=False, na=False)]

    # drop those rows where description=='description'
    scores = scores[scores['description'] != 'description']

    scores = scores[scores_types]
    scores = scores.apply(pd.to_numeric, errors='ignore')

    # Sort the DataFrame by the 'reweighted_sc' column and save it
    scores_sorted = scores.sort_values(by='reweighted_sc', ascending=True)

    # This is to reproduce previous behaviour and debugging
    scores_sorted.to_csv(output_file, sep='\t', index=False, header=False)

    return scores_sorted


def calculate_actual_radius(pose, radius=2):
    """
    Calculate radius to use based on the peptide length.
    
    :param pose: PyRosetta pose
    :param radius: The radius parameter as a float.

    :return: the calculated radius
    """
    
    # We assume two chains: receptor first, peptide second
    if len(get_chains(pose)) == 2:
        peptide_length =  len(pose.chain_sequence(2))
        total_length =  len(pose.chain_sequence(1)) + peptide_length
    else:
        # cluster only by peptide backbone
        peptide_length = len(pose.chain_sequence(1))

    # Calculate the actual radius using the formula
    actual_radius = math.sqrt(peptide_length / total_length) * radius

    # Output the actual radius
    print(f"Actual radius is {actual_radius}")

    return actual_radius


def run_clustering(all_poses, actual_radius, cluster_by_bb=True, max_num_clusters=None, outdir='clustering/'):
    """
    This runs the clustering (Phil style, the same as in the default cmdline application, and then sorts the structures according to the score defined in read_poses_from_silentfile.
    :param all_poses: list of PyRosetta Pose objects to be clustered
    :param actual_radius: calculated from peptide length
    :param cluster_by_bb: cluster by backbone or all-atom
    :param max_num_clusters: maximum number of clusters to create. default is the number of poses

    :return: the results of the clustering
    """
    clphil = ClusterPhilStyle()
    clphil.set_cluster_radius(actual_radius)

    #Set type of clustering
    if cluster_by_bb:
        clphil.set_cluster_by_protein_backbone(True)
    else:
        clphil.set_cluster_by_protein_backbone(True)

    for pose in all_poses:
        clphil.apply(pose)

    if not max_num_clusters:
        max_num_clusters = len(all_poses)
    clphil.do_clustering(max_num_clusters)
    clphil.do_redistribution()

    cl_list = clphil.get_cluster_list()
    clphil.sort_groups_by_energy()
    clphil.sort_each_group_by_energy()
    clphil.print_summary()
    clphil.print_cluster_assignment()

    # only if clustering/ directory exists
    if os.path.isdir('clustering/'):
        clphil.print_cluster_PDBs('clustering/')
    else:
        print("No clustering directory found, skipping printing cluster PDBs")

    return clphil


def copy_final_pdbs(cluster_id, member_id, rank):
    clustering_file = f'clustering/c.{cluster_id}.{member_id}.pdb'
    target_file = f'results/model_rank_{rank}.pdb'

    try:
        copyfile(clustering_file, target_file)
        print(f'Copied {clustering_file} to {target_file}')
    except:
        raise RuntimeError(f'Copying {clustering_file} to {target_file} was unsuccessful')


def process_clusters(clusters, scores=['I_sc', 'reweighted_sc', 'rmsBB_if']):
    """
    This function gets the top 10 lowest scoring clusters and outputs their cluster centers as solutions.
    Ranking according to Rosetta reweighted or other score, specified in read_poses_from_silentfile function.
    param: clusters from clustering

    output: files in clustering/ and results/ directories
    """
    # Gather values for the final file
    pose_list = []
    for index, cl in enumerate(clusters.get_cluster_list()):
        for p in range(cl.size()):
            pose = clusters.get_pose_list()[cl[p]]
            pose_scores = pose.scores
            pose_data = [tag_from_pose(pose), index, p]
            for score in scores:
                pose_data.append(round(pose_scores[score], 3))
                
            pose_list.append(pose_data)

    pose_df = pd.DataFrame(pose_list, columns=["Decoy_ID", "Cluster_no", "Member_ID"]+scores)
    pose_df = pose_df.apply(pd.to_numeric, errors='ignore')

    # Output dataframe similarly to previous runs for backward compatibility
    only_cluster_centers_df = pose_df.query('Member_ID == 0')

    return only_cluster_centers_df


def run_finalize(only_cluster_centers_df):
    """
    This will write out the dataframes and copy the top 10 structures to the results directory.
    """
    
    # Output dataframe similarly to previous runs for backward compatibility
    only_cluster_centers_df[["Decoy_ID", "Cluster_no", "Member_ID"]].to_csv('clustering/cluster_list',
                                                                            index=False, header=False, sep=' ')
    only_cluster_centers_df.head(10).to_csv(
        'clustering/cluster_list_reweighted_sc_sorted', index=False, header=False, sep=' ')
    only_cluster_centers_df.sort_values('I_sc').head(10).to_csv('clustering/cluster_list_I_sc_sorted',
                                                                index=False, header=False, sep=' ')
    
    # Create a df for the final 10 pdb-s
    top10_df = only_cluster_centers_df.sort_values('reweighted_sc').head(10)
    top10_df.insert(0, 'Rank', range(1, len(top10_df) + 1))
    top10_df.to_csv('results/final_scores', index=False, header=True, sep=' ')

    # Copy the structures to the results directory
    top10_df.apply(lambda row: copy_final_pdbs(row["Cluster_no"], row["Member_ID"], row["Rank"]), axis=1)


def clean_dir(dirs=[]):
    """
    Clean the directory from every file. Do not throw error if the file does not exist.
    """
    for dir in dirs:
        for file in os.listdir(dir):
            try:
                os.remove(file)
            except FileNotFoundError:
                pass

###########################

def get_elements_of_selection(sele):
	myspace = {'res_numbers': []}
	cmd.iterate(sele, 'res_numbers.append(resi)', space=myspace)  # iterate through model receptor interface
	output = np.unique(myspace['res_numbers'])
	
	return (output)


def find_shortest_chain(object_name):
	chain_lengths = {}
	
	# Iterate over all chains in the object and count residues
	chains = cmd.get_chains(object_name)
	for chain in chains:
		residues = get_elements_of_selection(f'chain {chain} and {object_name}')
		chain_lengths[chain] = len(set(residues))  # Use a set to count unique residues
	
	# Find the shortest chain
	shortest_chain = min(chain_lengths, key=chain_lengths.get)
	
	return chains[0], shortest_chain


def create_figures():
	pymol.pymol_argv = ['pymol', '-cq']
	cmd.loadall('model_rank*pdb')
	
	# load also native if exists
	native_pdbs = glob.glob('../*native*.pdb')
	if len(native_pdbs) > 0:
		native_pdb_file = native_pdbs[0]
		cmd.load(native_pdb_file, 'native')
		ch_rec_native, ch_pep_native = find_shortest_chain('native')
	else:
		native_pdb_file = None
		
	ch_rec_model, ch_pep_model = find_shortest_chain('model_rank_1')
	
	# Align first receptor in all models to native
	target = 'native' if native_pdb_file else 'model_rank_1'
	cmd.extra_fit(f"chain {ch_rec_model} and model*", target, "super")
	
	# First visualization: full complex
	cmd.hide("everything")

	cmd.show("cartoon", f"chain {ch_pep_model} and model*")
	cmd.show("cartoon", f"not chain {ch_pep_model} and model*")
	cmd.color("white", f"not chain {ch_pep_model} and model*")
	
	if native_pdb_file:
		ch_rec_native, ch_pep_native = find_shortest_chain('native')
		cmd.show("cartoon", f"chain {ch_pep_native} and native")
		cmd.spectrum("count", "rainbow_cycle", f"chain {ch_pep_native}")
		cmd.color("white", "native")
	
	cmd.bg_color("white")
	cmd.set("antialias", 2)
	cmd.set("ray_trace_mode", 1)
	cmd.orient()
	cmd.zoom()
	cmd.ray(1024)
	cmd.png("complex_view.png")
	
	# Subsequent visualizations for each model
	models = ["model_rank_1", "model_rank_2", "model_rank_3", "model_rank_4", "model_rank_5"]
	colors = ["cbag", "cbac", "cbay", "cbao", "cbap"]
	
	for model, color in zip(models, colors):
		try:
			cmd.hide("everything")
			cmd.show("cartoon", model)
			cmd.orient(f"chain {ch_pep_model}")
			cmd.show("sticks", f"chain {ch_pep_model} and {model} and not (name O or name C or name N)")
			cmd.hide("everything", f"chain {ch_pep_model} and {model} and hydrogen")
			getattr(cmd.util, color)(f"chain {ch_pep_model} and {model}")
			cmd.ray(1024)
			cmd.png(f"{model}.png")
		except:
			print(f'{model} does not exists')

#    cmd.save('models.pse')


############################################


def landscape_by_score(scores, score_type, rmsd_type):
	plt.style.use('seaborn-v0_8-whitegrid')
	f, ax = plt.subplots(1, 1, figsize=(10, 8))
	plt.scatter(scores[rmsd_type], scores[score_type], s=85, alpha=0.2)
	plt.ylabel(score_type, fontsize=20)
	plt.xlabel(rmsd_type, fontsize=20)
	plt.yticks(size=15)
	plt.xticks(size=15)
	plt.title(f"{rmsd_type} vs {score_type}", fontsize=25)
	ax.set_xlim(left=0)
	legend = plt.legend(frameon=3, loc="upper right", edgecolor='inherit', fontsize='20',
	                    title_fontsize='20')  # "upper left"
	frame = legend.get_frame()
	frame.set_color('white')
	f.savefig(f"results/{rmsd_type}_VS_{score_type}_landscape.png")


def save_score_file(scores_file, return_all_scores=True, score_type='reweighted_sc', rmsd_type='rmsBB_if'):
	if return_all_scores:
		out_columns = ['reweighted_sc', 'score', 'I_sc', 'pep_sc',
		               'pep_sc_noref', 'fa_atr', 'fa_rep', 'fa_sol',
		               'fa_dun', 'hbond_sc', 'rmsBB_if',
		               'rmsALL_if', 'rmsBB', 'startRMSbb', 'description']
	else:
		out_columns = [score_type, rmsd_type, 'description']
	
	df = scores_file[out_columns]
	df.to_csv("results/filtered_scores.tsv", sep="\t", index=False)
	
	return out_columns


def run_post_processing(scores, return_all_scores=True, score_type='reweighted_sc', rmsd_type='rmsBB_if'):
    n_rows = int(len(scores) * 0.9)
    
    # sorted and extract pdb files by score, return only 90% of the rows
    scores_90 = scores.sort_values(score_type, ascending=True).head(n_rows)  
    landscape_by_score(scores_90, score_type, rmsd_type)
    save_score_file(scores_90, return_all_scores, score_type, rmsd_type)  # save results in tsv file
	
    os.chdir('results')
    create_figures()


            
def parse_args():
    parser = argparse.ArgumentParser(description="Run clustering and finalization steps.")
    parser.add_argument("-w", "--work_dir", help="Directory containing decoys and other files.", default=".")
    parser.add_argument("-s", "--score_file", help="Path to the score", default="score.sc")
    parser.add_argument("-d", "--silent_file", help="Path to the score", default="decoys.silent")
    parser.add_argument("-t", "--tags_to_remove", nargs="+", default=[], help="Tags to exclude from clustering.")
    parser.add_argument("-r", "--radius", type=float, default=2.0, help="RMSD threshold for clustering.")
    parser.add_argument("-x", "--topX", type=int, default=0.01, help="Ratio of decoys to take for clustering (0<x<=1) or exact number x>1")
    parser.add_argument("--no-rescore", action="store_true", help="Do not rescore against top scoring structure.")

    return parser.parse_args()


def main():
    pyrosetta.init('-mute FlexPepDockingProtocol FlexPepDockingPoseMetrics protocols.flexPepDocking core.pack core.scoring basic.io.database')

    args = parse_args()

    if not os.path.isdir(args.work_dir):
        raise RuntimeError('Working directory does not exist')
    else:
        os.chdir(args.work_dir)

    # Make clustering & results directory, dont throw error if it exists
    os.makedirs("clustering", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    clean_dir(['clustering', 'results'])

    if not os.path.isfile(args.silent_file):
        raise RuntimeError(f'Silent file {args.silent_file} does not exist.')

    if not os.path.isfile(args.score_file):
        raise RuntimeError(f'Score file {args.score_file} does not exist.')

    # First sort the scores and select and load top 1% (or other specified)
    scores_sorted = extract_and_sort_score_file(score_file=args.score_file, output_file='sorted.sc')
    tags = filter_and_select_top_structures(scores_sorted, tags_to_remove=args.tags_to_remove, topX=args.topX).description.tolist()

    if not args.no_rescore:
        scores_sorted, all_poses = rescore_models(args.silent_file, scores_sorted, tags, output_file='rescore.sc')
        scores_sorted.to_csv('sorted.sc', sep='\t', index=False, header=False)
    else:
        all_poses = read_poses_from_silentfile(args.silent_file, tags)
    actual_radius = calculate_actual_radius(all_poses[0], args.radius)

    # Run clustering and process the results
    clusters = run_clustering(all_poses, actual_radius)
    only_cluster_centers_df = process_clusters(clusters)
    run_finalize(only_cluster_centers_df)
    
    # Run post-processing to plot energy landscapes and save cleaned score file
    run_post_processing(scores_sorted, return_all_scores=True, score_type='reweighted_sc', rmsd_type='rmsBB_if')


def paired_residue_inds(a, b):
    """Get paired indices of common residues in two structures."""
    # From: https://gist.github.com/asford/c2404c8b045700f016fda8893325c807 by Alex Ford
    aseq = packed_pose.to_pose(a).sequence()
    bseq = packed_pose.to_pose(b).sequence()

    astart, bstart, align_len = difflib.SequenceMatcher(a=aseq, b=bseq).find_longest_match(0, len(aseq), 0, len(bseq))

    return [(astart + i, bstart + i) for i in range(align_len)]


def aligned_decoys(decoy_collection, align_to=None):
    """CA-align all structures in collection to target (or first) structure."""
    # From: https://gist.github.com/asford/c2404c8b045700f016fda8893325c807 by Alex Ford

    decoy_collection = packed_pose.to_dict(decoy_collection)

    if align_to is None:
        align_to = decoy_collection[0]
    align_to_pose = packed_pose.to_pose(align_to)

    aligned_decoys = []
    for decoy in decoy_collection:
        work_pose = packed_pose.to_pose(decoy)
        ca_map = map_core_id_AtomID_core_id_AtomID()
        for w, a in paired_residue_inds(work_pose, align_to_pose):
            ca_map[AtomID(work_pose.residue(w+1).atom_index("CA"), w+1)] = AtomID(
                align_to_pose.residue(a+1).atom_index("CA"), a+1
            )

        superimpose_pose(work_pose, align_to_pose, ca_map)
        aligned_decoys.append(work_pose)

    return aligned_decoys

if __name__ == "__main__":
    main()
