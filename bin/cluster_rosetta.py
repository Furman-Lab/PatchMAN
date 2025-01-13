

import os
import pandas as pd
import sys
import subprocess
import math
import os


import pyrosetta
# Initialize PyRosetta
from pyrosetta import *
from pyrosetta.rosetta.protocols.multistage_rosetta_scripts.cluster import ClusterMetric
from pyrosetta.rosetta.protocols.cluster import Cluster
from pyrosetta.rosetta.core.io.silent import SilentFileData, SilentStruct,SilentFileOptions
from pyrosetta.io.silent_file_map import SilentFileMap
pyrosetta.init()


# Initialize PyRosetta

def cluster_structures(
    poses=None,
    silent_file=None,
    cluster_radius=2.0,
    fullatom=True,
    tags=None
):
    """
    Clusters structures based on RMSD, supports fullatom and tags filtering.

    Args:
        poses (list): List of PyRosetta Pose objects.
        silent_file (str): Path to the silent file (optional).
        cluster_radius (float): RMSD threshold for clustering.
        fullatom (bool): If True, process silent file in full-atom mode.
        tags (list): Specific tags from the silent file to include in clustering.

    Returns:
        list: A list of clusters, each containing indices of the clustered poses.
    """
    # Initialize Cluster
    cluster = Cluster()
    #cluster.set_radius(cluster_radius)

    # Create the ClusterMetric object for scoring the clusters
    cluster_metric = ClusterMetric()
    all_poses = []

    # Add poses from the input list
    if poses:
        for pose in poses:
            cluster_metric.add_pose(pose)
        all_poses.extend(poses)

    # Process silent file if provided
    if silent_file:
        # Create SilentFileOptions and set fullatom flag
        options = SilentFileOptions()
        options.in_fullatom(fullatom)  # Ensure full-atom mode is respected

        silent_data = SilentFileData(silent_file, options)
        silent_file_map = SilentFileMap(silent_file)
        print(silent_data)

        # Loop through each structure in the silent file and extract its tag
        for tag, silent_struct in silent_file_map.items():
            print(f"Processing tag: {tag}")

            # Access the Pose object from the SilentStructure

            if tag not in tags:
                continue  # Skip if tag is not i
            #silent_struct = silent_file_map[tag]  # Retrieve the SilentStructure for the given tag
            # Initialize a Pose object

            #silent_struct = silent_data.get_structure(tag)
            pose = Pose()
            silent_struct.fill_pose(pose)
         #   cluster.add_pose(pose)
            all_poses.append(pose)

    if not all_poses:
        raise ValueError("No structures provided. Supply either poses or a silent file.")
    score_function = get_score_function()

    # Perform clustering
    cluster.cluster_poses(all_poses, cluster_radius, score_function)

    clusters = cluster_metric.cluster(cluster_radius)

    # Report clusters
    cluster_results = []
    for i, cluster in enumerate(clusters):
        print(f"Cluster {i}: {len(cluster)} structures")
        cluster_results.append(cluster)

    return cluster_results, all_poses



def extract_and_sort_score_file(score_file, output_file="sorted.csv"):
    path = os.path.dirname(score_file)
    if not '/' in output_file:
        output_file=path+'/'+output_file
    """
    Extract columns from the score file (score.sc), sort by 'reweighted_sc', and return the result as a pandas DataFrame.

    :param score_file: Path to the score file (e.g., "score.sc").
    :param output_file: Path to the output file (e.g., "sorted.sc"). This file will be created only if not already exists.
    :return: pandas DataFrame containing the sorted data.
    """
    # If output_file exists, return the corresponding DataFrame
    if not os.path.isfile(output_file):
         df= pd.read_csv(output_file, delim_whitespace=True)
         df[[ 'I_sc','reweighted_sc', 'rmsBB', 'rmsBB_if']] = df[['I_sc','reweighted_sc', 'rmsBB', 'rmsBB_if']].astype(float)
         df[ 'reweighted_sc'] = df['reweighted_sc'].astype(float)
         df = df.sort_values(by='reweighted_sc', ascending=True)
         df = df.drop('SCORE:', axis=1)

         return df

    # List of scores to extract
    scores = ['reweighted_sc', 'I_sc', 'rmsBB', 'rmsBB_if', 'description']
    """
    Processes the Rosetta score file and extracts the specified scores into a pandas DataFrame.

    :param score_file: Path to the Rosetta score file.
    :param scores: List of score names to extract (e.g., ['reweighted_sc', 'I_sc']).
    :return: pandas DataFrame containing the extracted scores.
    """
    # Initialize variables for processing the file
    header = []

    # Open and process the score file
    with open(score_file, 'r') as file:
            # Read the first line for the header
            p=1
            while p<3 :
                   header = file.readline().strip().split('\t')
                   p=p+1
                   if 'total_score' in header[0]:
                        header=header[0].split(' ')
                        header = [item for item in header if item != '']  # Remove empty strings
                        break

    # Now read the rest of the file into a DataFrame
    df = pd.read_csv(score_file, sep='\s+', header=None, names=header, skiprows=1)
    df = df[~df['SCORE:'].str.contains('SEQUENCE:', case=False, na=False)]
    df = df.drop('SCORE:', axis=1)
    # Create a DataFrame with the extracted data
    df[['I_sc', 'reweighted_sc', 'rmsBB', 'rmsBB_if']] = df[['I_sc', 'reweighted_sc', 'rmsBB', 'rmsBB_if']].astype(float)
    df['reweighted_sc'] = df['reweighted_sc'].astype(float)

    # Sort the DataFrame by the 'reweighted_sc' column
    df_sorted = df.sort_values(by='reweighted_sc', ascending=True)
    # If output_file doesn't exist, save the DataFrame to the file
    if not os.path.isfile(output_file):
        df_sorted.to_csv(output_file, sep=' ', index=False, header=True)

    return df_sorted


def filter_and_select_top_structures(df,discription_remove=[], percent=1):
    """
    Filters the score file to exclude certain descriptions and selects top X percent structures based on score.

    :param score_file: Path to the score file (e.g., "sorted.sc").
    :param percent: Percentage of top structures to select (default is 1%).
    :return: A DataFrame containing the filtered top X percent structures.
    """
    # Read the score file into a DataFrame

    # Assuming columns are well defined, for example: ['score', ..., 'description']
    # Filter possible disulfide bridges: process based on 'description' and 'score'

    descriptions_to_remove = []
    for _, row in df.iterrows():
        # Extract the score for each row
        score = row['reweighted_sc']
        description = row['description']
        # If the score is less than -1000000, mark for removal
        if score < -1000000:
            descriptions_to_remove.append(description)
        else:
            break  # Stop processing further once a score >= -1000000 is encountered

    descriptions_to_remove = descriptions_to_remove+discription_remove
    # Filter out rows that match the descriptions to be removed
    if descriptions_to_remove:
        # Create a regex pattern from descriptions to remove
        all_filter = '|'.join(descriptions_to_remove)
        df_filtered = df[~df['description'].str.contains(all_filter, regex=True, na=False)]
    else:
        df_filtered = df.copy()

    # Calculate the number of rows to take for the top X percent
    num_lines = len(df_filtered)
    if num_lines > 1000:
        num_str = num_lines // 100
    elif num_lines > 250:
        num_str = num_lines // 10
    else:
        num_str = num_lines // 2

    # Select the top X percent based on sorted scores
    top_percent_df = df_filtered.nsmallest(num_str, 'reweighted_sc')
    top_percent_df.to_csv('top_percent_df.csv')

    return top_percent_df


def actual_radius(peptide_length, rec_len, radius):
    """
    Run clustering based on peptide sequence and radius.

    :param peptide_sequence: The peptide sequence as a string.
    :param radius: The radius parameter as a float.
    """
    # Create a new directory for clustering
    clustering_dir = "clustering"
    os.makedirs(clustering_dir, exist_ok=True)

    # Find the first PDB file that matches the pattern
    pdb_files = [f for f in os.listdir('.') if f.endswith('_0001.pdb')]
    pdb_file = pdb_files[0] if pdb_files else None

    if pdb_file is None:
        print("No PDB files found.")
        return

    # Calculate the actual radius using the formula
    actual_radius = math.sqrt(peptide_length / rec_len) * radius
    # Output the actual radius
    print(f"Actual radius is {actual_radius}")

    return actual_radius


def run_clustering_with_function(dirr, score_file, peptide_length, rec_len, removed_tags=[], radius=2.0):

    """
    Run clustering using the `cluster_structures_with_rmsd` function.

    Parameters:
        dirr (str): Directory containing decoys and other files.
        score_file (str): Path to the score file.
        peptide_length (int): Length of the peptide.
        rec_len (int): Length of the receptor.
        removed_tags (list): Tags to exclude from clustering.
        radius (float): RMSD threshold for clustering.

    Returns:
        tuple: Clustering results (clusters, top_percent_df).
    """
    # Change working directory
    os.chdir(dirr)

    # Get the sorted file
    sorted_file = extract_and_sort_score_file(score_file)

    # Filter top percent structures
    top_percent_df = filter_and_select_top_structures(sorted_file, discription_remove=removed_tags)
    top_percent_descriptions = top_percent_df['description'].tolist()

    # Calculate actual radius
    actualR = actual_radius(peptide_length, rec_len, radius)

    # Silent file containing decoys
    decoys_file = dirr+ "/decoys.silent"

    # Use cluster_structures_with_rmsd() for clustering
    print(top_percent_descriptions)
    clusters = cluster_structures(
        silent_file=decoys_file,
        tags=top_percent_descriptions,
        cluster_radius=actualR
    )#pyrossta code for custring

    # Return clusters and top_percent_df
    return clusters, top_percent_df



def run_clustering(dirr, score_file, peptide_length, rec_len, removed_tags=[],redius=2):
    # File paths
    os.chdir(dirr)
    ROSETTA_BIN = "/vol/ek/share/rosetta/rosetta_src_2020.28.61328_bundle/rosetta.source.release-260/main/source/bin/cluster.default.linuxgccrelease"
    decoys_file = os.path.join(dirr, "decoys.silent")  # Assuming decoys.silent is in the provided directory
    # top1percent_file = os.path.join(dirr, "top1percent")

    # Get the sorted file
    sorted_file = extract_and_sort_score_file(score_file)
    # Filter top percent structures
    top_percent_df = filter_and_select_top_structures(sorted_file, discription_remove=removed_tags)
    top_percent_descriptions = top_percent_df['description'].tolist()
    # Convert descriptions to a string with space-separated values (this is passed to the Rosetta command)
    descriptions_str = " ".join(top_percent_descriptions)
    # Calculate actual radius
    actualR = actual_radius(peptide_length, rec_len, redius)

    # Prepare the sbatch command
    sbatch_command = f"""
    #!/bin/bash
    #SBATCH --time=93:00:00
    #SBATCH --mem=2000m
    #SBATCH -job-name=clustering 
	#SBATCH	--nice=8000 
	#SBATCH	--kill-on-invalid-dep=yes

    {ROSETTA_BIN} -in:file:silent {decoys_file} -in:file:silent_struct_type binary -cluster:radius {actualR} -in:file:fullatom -tags {descriptions_str} > clog
    """
    # Create a temporary bash script to submit
    sbatch_script = os.path.join(dirr, "cluster_sbatch.sh")
    with open(sbatch_script, 'w') as f:
        f.write(sbatch_command)

    # Run sbatch for clustering
    try:
        clustering_jid = subprocess.check_output(["sbatch", sbatch_script]).decode('utf-8').split()[-1]
        print(f"Clustering job submitted successfully with job ID {clustering_jid}")
    except subprocess.CalledProcessError as e:
        print(f"Error while submitting clustering job: {e}")
        return None
    # Optionally remove the temporary sbatch script
    # os.remove(sbatch_script)
    # Return clustering job ID for further dependency handling
    return clustering_jid, top_percent_df


def process_clog_and_top1percent(clog_file, top1percent_df):
    # Get the number of rows from the top1percent DataFrame (equivalent to `wc -l` in bash)
    x = len(top1percent_df)
    # Read the last `x + 5` lines from the `clog` file (equivalent to `tail -$((x+5)) clog`)
    with open(clog_file, 'r') as f:
        clog_lines = f.readlines()
    # Extract lines from -x to -5
    last_x_lines = clog_lines[-(x + 5):-5]

    # Filter lines containing '_0001' (equivalent to `grep _0001`)
    filtered_lines = [line for line in last_x_lines if '_0001' in line]

    # Process each line: extract the matching columns (equivalent to `awk` and `print $i,$(i+1),$(i+2)`)
    processed_data = []
    for line in filtered_lines:
        columns = line.split()
        for i in range(len(columns)):
            if '_0001' in columns[i]:
                # Extract columns i, i+1, i+2
                processed_data.append([columns[i], columns[i + 1], columns[i + 2]])

    # Convert the processed data into a DataFrame
    df = pd.DataFrame(processed_data, columns=['description', 'Cluster_no', 'Member_ID'])

    # Remove any empty rows (equivalent to `sed '/^$/d'`)
    df.dropna(inplace=True)

    # Output the DataFrame
    return df


def process_cluster_and_score(cluster_list, df_score, pdb_columns='description', output_dir="results"):
    """
    Process the cluster list DataFrame and score DataFrame to extract and combine relevant information.

    Parameters:
    - cluster_list (pd.DataFrame): DataFrame containing cluster list information.
    - df_score (pd.DataFrame): DataFrame containing score information (score.sc file).
    - pdb_columns (str): The column name used for matching between cluster_list and df_score (default is 'description').
    - output_dir (str): Directory to save the final processed results. Default is "results".

    Returns:
    - pd.DataFrame: A DataFrame containing the extracted results.
    """
    # Create an empty list to store the results
    results = []

    # Loop over each value in the specified column of cluster_list
    for i in cluster_list[pdb_columns]:
        # Find the matching rows in df_score based on the value `i` (equivalent to `grep $i` in Bash)
        matching_rows = df_score[df_score[pdb_columns] == i]

        # Assuming you want to add the first row of the matching rows (equivalent to `head -1` in Bash)
        if not matching_rows.empty:
            first_row = matching_rows.iloc[0]  # Get the first match

            # Extract the relevant columns (I_sc, reweighted_sc, rmsBB_if)
            relevant_values = first_row[['I_sc', 'reweighted_sc', 'rmsBB_if']].values.tolist()

            # Append the extracted values to the results list
            results.append(relevant_values)

    # Convert the results into a DataFrame
    results_df = pd.DataFrame(results, columns=['I_sc', 'reweighted_sc', 'rmsBB_if'])

    # Combine with the original cluster_list DataFrame
    cluster_list_combined = pd.concat([cluster_list, results_df], axis=1)

    # Save the final result to the specified output directory
    cluster_list_combined.to_csv(f'{output_dir}/cluster_list_sc.csv', index=False)

    print(f"Results saved to {output_dir}/cluster_list_sc.csv")

    return cluster_list_combined


def sort_and_combine_clusters(cluster_list, pdb_list_sc, final_results_path):
    """
    Sort and combine the cluster list DataFrame with the score DataFrame.
    The combined results will be sorted by `I_sc` and `reweighted_sc`,
    keeping only the top 10 entries for each.

    Parameters:
    - cluster_list (pd.DataFrame): DataFrame containing cluster information.
    - pdb_list_sc (pd.DataFrame): DataFrame containing score information (I_sc, reweighted_sc, rmsBB_if).
    - final_results_path (str): Directory path to store the final results.

    Returns:
    - pd.DataFrame: Combined and sorted DataFrame based on I_sc and reweighted_sc.
    """
    # Combine cluster_list and pdb_list_sc into one DataFrame
    cluster_list_combined = pd.concat([cluster_list, pdb_list_sc], axis=1)

    # Sort by I_sc column and keep unique Cluster_no, selecting top 10 entries
    cluster_list_I_sc_sorted = cluster_list_combined.sort_values(by=['I_sc'], ascending=True)
    cluster_list_I_sc_sorted = cluster_list_I_sc_sorted.drop_duplicates(subset=['Cluster_no'], keep='first').head(10)

    # Sort by reweighted_sc column and keep unique Cluster_no, selecting top 10 entries
    cluster_list_reweighted_sc_sorted = cluster_list_combined.sort_values(by=['reweighted_sc'], ascending=True)
    cluster_list_reweighted_sc_sorted = cluster_list_reweighted_sc_sorted.drop_duplicates(subset=['Cluster_no'],
                                                                                          keep='first').head(10)

    # Prepare the final sorted DataFrames
    cluster_list_I_sc_sorted = cluster_list_I_sc_sorted[
        ['Decoy_ID', 'Cluster_no', 'Member_ID', 'I_sc', 'reweighted_sc', 'rmsBB_if']]
    cluster_list_reweighted_sc_sorted = cluster_list_reweighted_sc_sorted[
        ['Decoy_ID', 'Cluster_no', 'Member_ID', 'I_sc', 'reweighted_sc', 'rmsBB_if']]

    # Creating output directories
    os.makedirs(final_results_path, exist_ok=True)

    # Save the sorted results to CSV (or any other format)
    cluster_list_I_sc_sorted.to_csv(f"{final_results_path}/cluster_list_I_sc_sorted.csv", index=False)
    cluster_list_reweighted_sc_sorted.to_csv(f"{final_results_path}/cluster_list_reweighted_sc_sorted.csv", index=False)
    cluster_list_I_sc_sorted.to_csv(f"{final_results_path}/final_scores.csv", index=False)
    cluster_list_reweighted_sc_sorted.to_csv(f"{final_results_path}/results_info.csv", index=False)
    # Copy results to final results directory


    # Mimicking the copying of .pdb files for top 10 entries
    for _, row in cluster_list_reweighted_sc_sorted.iterrows():
        pdb_filename = f"c.{row['Cluster_no']}.{row['Member_ID']}.pdb"
        os.system(f"cp clustering/{pdb_filename} {final_results_path}")#Dont think its good to do it like this

    return cluster_list_I_sc_sorted, cluster_list_reweighted_sc_sorted


def run_finalize(clusters,top1percent_df, dirr, score_file,output_file="sorted.sc"):
    cluster_list=process_clog_and_top1percent(clusters, top1percent_df)
    sorted_scores_df=extract_and_sort_score_file(score_file,output_file )
    cluster_list_combined=process_cluster_and_score(cluster_list,sorted_scores_df)
    final_cluster_sorted, final_reweighted_sorted = sort_and_combine_clusters(cluster_list_combined, sorted_scores_df,
                                                                              dirr)
    return final_cluster_sorted, final_reweighted_sorte

def parse_args(args):
    parser = argparse.ArgumentParser(description="Run clustering and finalization steps.")
    parser.add_argument("--dir", required=True, help="Directory containing decoys and other files.")
    parser.add_argument("--score_file", required=True, help="Path to the score file.")
    parser.add_argument("--peptide_length", type=int, required=True, help="Length of the peptide.")
    parser.add_argument("--removed_tags", nargs="+", default=[], help="Tags to exclude from clustering.")
    parser.add_argument("--radius", type=float, default=2.0, help="RMSD threshold for clustering.")
    return parser.parse_args(args)


def main():
    #only my example
    dirr = "/sci/nosnap/fora/nirit/AF2_V3_LNR/run_Pipline_iptmCutOff_Focus_PatchMan/results/ALL/recycel9/1awr/AF2_pdbs_FPD/"
    score_file = "/sci/nosnap/fora/nirit/AF2_V3_LNR/run_Pipline_iptmCutOff_Focus_PatchMan/results/ALL/recycel9/1awr/AF2_pdbs_FPD/score.sc"
    peptide_length = 10
    removed_tags = []#["tag1", "tag2"]
    radius = 2.0

    clusters, top_percent_df = run_clustering_with_function(
        dirr, score_file, peptide_length, rec_len, removed_tags, radius
    )
    run_finalize(clusters, top_percent_df, dirr, score_file, output_file="sorted.sc")
    print(clusters,top_percent_df)
    return clusters,top_percent_df


clusters,top_percent_df=main()

print(clusters,top_percent_df)
