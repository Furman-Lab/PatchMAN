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


def create_tasks(list_of_inputs, dict_of_options, nstruct=1):
    for input in list_of_inputs:
        basename = os.path.splitext(os.path.basename(input))[0]
        for index in range(1, int(nstruct) + 1):
            yield {
            "options": "-ex1",
            "extra_options": dict_of_options,
            "set_logging_handler": "interactive",
            "input": input,
            "numbered_basename": f"{basename}_{index:04d}",
            "keep_best": True           # this is okay for runs with equal length of peptide
            }


def protocol(pose, **kwargs):
    pose = io.pose_from_file(kwargs["input"]) # for some reason, the pose is not passed correctly
    xml = f"""
        <ROSETTASCRIPTS>
            <SCOREFXNS>
            </SCOREFXNS>
            <MOVERS>
                <FlexPepDock name="fpd" pep_refine="1" extra_scoring="1" lowres_preoptimize="1" {kwargs["fpd_args"]}/>
            </MOVERS>

            <PROTOCOLS>
                <Add mover_name="fpd" />
            </PROTOCOLS>
        </ROSETTASCRIPTS>
    """

    # Name the decoys according to the usual Rosetta naming scheme. Need to convert to regular Pose for this.
    packed_pose = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)(pose.pose.clone())
    out_pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
    out_pose.pdb_info().name(kwargs['numbered_basename'])

    # Add extra scores to Pose object
    extra_scores = dict(get_string_real_pairs_from_current_job())
    for k, v in extra_scores.items():
        if not k.startswith('best') and not kwargs['keep_best']:
            out_pose.scores[k] = v

    # Save into the silent file
    io.to_silent(out_pose, f"{kwargs['PyRosettaCluster_output_path']}/output.silent")

    return


def run_fpd_cluster(list_of_inputs, dict_of_options, output_path='.', nstruct=1):
    client = Client()

    PyRosettaCluster(
        tasks=create_tasks(list_of_inputs, dict_of_options, nstruct),
        client=client,
        scratch_dir=output_path,
        output_path=output_path,
        nstruct=1,                      # here it is only one, as otherwise the decoy numbering cannot be tracked
        sha1=None,                      # to prevent an error with git
        scorefile_name="scores.json",   # with dry-run, this wont be generated anyway
        dry_run=True,
        decoy_dir_name='.',             # dont create  unnecessary directories
        logs_dir_name='.'               # dont create unnecessary directories
    ).distribute(protocols=[protocol])

    # Make a score file: open the silent file and extract lines starting with 'SCORE'
    with open(os.path.join(output_path, "decoys.silent"), "r") as infile:
        score_lines = [line for line in infile if line.startswith("SCORE")]

    # Write the extracted lines to a new file
    with open(os.path.join(output_path, "scores.sc"), "w") as outfile:
        outfile.writelines(score_lines)


def main():
    parser = argparse.ArgumentParser(description="Controls the FlexPepDock refinement part of the protocol.")
    parser.add_argument("-t", "--nstruct", type=int, help="Number of FlexPepDock decoys to generate. Default: 1")
    parser.add_argument("-m", "--min_rec_bb", action="store_true", help="Use receptor backbone minimization? Doubles runtime. Default: False")
    parser.add_argument("-u", "--unboundrot", help="Add unbound rotamers for the receptor. Default: None")
    parser.add_argument("-a", "--native", help="Native structure for comparison. Default: None")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Enable benchmarking mode")
    parser.add_argument("-r", "--redundancy", action="store_true", help="Filter for template redundancy. Default: False")
    parser.add_argument("-c", "--cpu", type=int, help="Number of CPUs to use. Default: All available")
    args = parser.parse_args()

    # Get number of threads to use
    if args.cpu:
        procs = args.cpu
    elif "SLURM_NPROCS" not in os.environ:
        procs = os.cpu_count()
    else:
        procs = int(os.environ["SLURM_NPROCS"])

    print(f"Running FlexPepDock refinement on {procs} threads")

    # For init Rosetta
    dict_of_init_opts = {
        f"-multithreading:total_threads": procs,
        f"-overwrite": "",
    }

    # Parse arguments
    if args.unboundrot:
        dict_of_init_opts["-unboundrot"] = args.unboundrot
    if args.native:
        dict_of_init_opts["-native"] = args.native

    fpd_args = ""
    if args.min_rec_bb:
        fpd_args += ' min_receptor_bb="1"'


    # Collect input files
    input_files = glob.glob("*0001.pdb")

    # Filter for peptide backbone templates
    if args.redundancy:
        input_files = keep_only_redundant_templates(input_files)

    # Benchmarking mode
    if args.benchmark:
        # Read the benchmark file and input list
        benchmark_names = Path("benchmark_file").read_text().splitlines()

        # Filter out elements from the input list that contain any benchmark name
        input_files = [line for line in input_files if
                         not any(name.lower() in line.lower() for name in benchmark_names)]

    # Default parameters
    n_templates = len(input_files)
    min_nstruct = n_templates * 10
    max_nstruct = 20000

    if args.nstruct is None:
        print("Nstruct needs to be set automatically")
        if min_nstruct < max_nstruct:
            nstruct = 10
        elif n_templates > max_nstruct:
            nstruct = 1
        else:
            nstruct = max_nstruct // n_templates
    else:
        nstruct = args.nstruct

    print(f"Refinement will generate {nstruct} decoys per input structure")

    # Clean up existing output files before running FlexPepDock
    for file in ["score.sc", "decoys.silent"]:
        Path(file).unlink(missing_ok=True)

    # Run FlexPepDock
    run_fpd_cluster(input_files, dict_of_init_opts, os.getcwd(), int(nstruct))

if __name__ == "__main__":
    main()


