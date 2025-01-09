import argparse
import glob
import os
import subprocess
from pathlib import Path
from collections import Counter

import pyrosetta
import pyrosetta.distributed.io as io
import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
from dask.distributed import Client
from pyrosetta import *
from pyrosetta.distributed.cluster import PyRosettaCluster
from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job


# pip3 install distributed tools dask attrs billiard

def keep_only_redundant_templates(input_files):

    # Extract the middle part (2nd to 4th segments) from filenames
    segments = []
    for filename in input_files:
        parts = filename.split('_')
        if len(parts) > 3:
            segments.append('_'.join(parts[1:4]))

    # Count occurrences of each segment
    counts = Counter(segments)

    # Get segments that occur more than once
    results = [key for key, count in counts.items() if count > 1]

    return results


def create_tasks(list_of_inputs, dict_of_options, nstruct, min_rec_bb):
    for input in list_of_inputs:
        basename = os.path.splitext(os.path.basename(input))[0]
        for index in range(1, nstruct + 1):
            print(input, basename, "dict", dict_of_options, f"{basename}_{index:04d}")
            yield {
            "options": "-ex1 -ex2aro -use_input_sc -overwrite",
            "extra_options": dict_of_options,
            "set_logging_handler": "interactive",
            "s": input,
            "numbered_basename": f"{basename}_{index:04d}",
            "min_rec_bb": min_rec_bb
            }


def protocol(pose, **kwargs):
    import pyrosetta
    import pyrosetta.distributed.io as io
    import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job

    pose = io.pose_from_file(kwargs["s"])

    xml = f"""
        <ROSETTASCRIPTS>
            <SCOREFXNS>
            </SCOREFXNS>
            <MOVERS>
                <FlexPepDock name="fpd" pep_refine="1" extra_scoring="1" lowres_preoptimize="1" min_receptor_bb='{str(int(kwargs["min_rec_bb"]))}'/>
            </MOVERS>

            <PROTOCOLS>
                <Add mover_name="fpd" />
            </PROTOCOLS>
        </ROSETTASCRIPTS>
    """

    # Name the decoys according to the usual Rosetta naming scheme. Need to convert to regular Pose
    packed_pose = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)(pose.pose.clone())
    out_pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
    out_pose.pdb_info().name(kwargs['numbered_basename'])

    # Add extra scores to Pose object
    extra_scores = dict(get_string_real_pairs_from_current_job())
    for k, v in extra_scores.items():
        if not k.startswith('best'):
            out_pose.scores[k] = v

    # Save into the silent file
    io.to_silent(out_pose, f"{kwargs['PyRosettaCluster_output_path']}/decoys.silent")

    return out_pose  # otherwise score file is not written


def run_fpd_cluster(client, list_of_inputs, dict_of_options, min_rec_bb=False, output_path='.', nstruct=1):
    
    PyRosettaCluster(
        tasks=create_tasks(list_of_inputs, dict_of_options, nstruct, int(min_rec_bb)),
        client=client,
        scratch_dir=output_path,
        output_path=output_path,
        nstruct=1,  # here it is only one, as otherwise the decoy numbering cannot be tracked
        sha1=None,  # to prevent an error with git
        scorefile_name="scores.json",  # with dryrun, this wont be generated
        dry_run=True,
        decoy_dir_name='.',  # dont create  unnecessary directories
        logs_dir_name='.'  # dont create unnecessary directories
    ).distribute(protocols=[protocol])

    # Make a score file: open the silent file and extract lines starting with 'SCORE'
    with open(os.path.join(output_path, "decoys.silent"), "r") as infile:
        score_lines = [line for line in infile if line.startswith("SCORE")]

    # Write the extracted lines to a new file
    with open(os.path.join(output_path, "scores.sc"), "w") as outfile:
        outfile.writelines(score_lines)


def set_nstruct(n_templates, min_nstruct, max_nstruct=20000, args=None):
    if args and args.nstruct:
        nstruct = args.nstruct
    else:
        if min_nstruct * n_templates < max_nstruct:
            nstruct = 10
        elif n_templates > max_nstruct:
            nstruct = 1
        else:
            nstruct = max_nstruct // n_templates

    print(f"Refinement will generate {nstruct} decoys per input structure for {n_templates} inputs")

    return nstruct


# Calculate the number of processors to use
def set_procs(args, nstruct, input_files):
    if args.cpu:
        procs = args.cpu
    elif "SLURM_NPROCS" not in os.environ:
        if nstruct * len(input_files) < os.cpu_count():
            procs = nstruct * len(input_files)
        else:
            procs = os.cpu_count()
    else:
        procs = int(os.environ["SLURM_NPROCS"])

    print(f"Running FlexPepDock refinement on {procs} threads")

    return procs


def clean_before_run():
    for file in ["score.sc", "decoys.silent"]:
        Path(file).unlink(missing_ok=True)


# process redundancy and/or benchmarking
def filter(input_files, redundancy, benchmark):
    if redundancy:
        input_files = keep_only_redundant_templates(input_files)

    if benchmark:
        # Read the benchmark file and input list
        benchmark_names = Path("benchmark_file").read_text().splitlines()

        # Filter out elements from the input list that contain any benchmark name
        input_files = [line for line in input_files if
                         not any(name.lower() in line.lower() for name in benchmark_names)]

    return input_files


# create a dictionary of options to init Rosetta. These need to be set from the cmdline
def create_init_opts_dict(procs, args):
    dict_of_init_opts = {
        f"-multithreading:total_threads": str(procs),
        f"-overwrite": "",
    }

    # Add arguments to the dict
    if args.unboundrot:
        dict_of_init_opts["-unboundrot"] = args.unboundrot
    if args.native:
        dict_of_init_opts["-native"] = args.native

    return dict_of_init_opts


def main():
    parser = argparse.ArgumentParser(description="Controls the FlexPepDock refinement part of the protocol.")
    parser.add_argument("-t", "--nstruct", type=int, help="Number of FlexPepDock decoys to generate. Default: 1")
    parser.add_argument("-m", "--min_rec_bb", action="store_true", help="Use receptor backbone minimization? Doubles runtime. Default: False")
    parser.add_argument("-u", "--unboundrot", type=str, help="Add unbound rotamers for the receptor. Default: None")
    parser.add_argument("-a", "--native", type=str, help="Native structure for comparison. Default: None")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Enable benchmarking mode")
    parser.add_argument("-r", "--redundancy", action="store_true", help="Filter for template redundancy. Default: False")
    parser.add_argument("-c", "--cpu", type=int, help="Number of CPUs to use. Default: All available")
    args = parser.parse_args()

    # Collect input files
    input_files = glob.glob("*0001.pdb")

    # Filter input files
    input_files = filter(input_files, args.redundancy, args.benchmark)

    # Set number of decoys to generate
    nstruct = set_nstruct(len(input_files), 10, 20000, args)

    # Set number of processors to use
    procs = set_procs(args, nstruct, input_files)

    # For init Rosetta
    dict_of_init_opts = create_init_opts_dict(procs, args)

    # Clean up existing output files before running FlexPepDock
    clean_before_run()

    # Run FlexPepDock
    client = Client()
    run_fpd_cluster(client, input_files, dict_of_init_opts, args.min_rec_bb, os.getcwd(), int(nstruct))

if __name__ == "__main__":
    main()


