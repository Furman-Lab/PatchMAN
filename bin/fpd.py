import argparse
import glob
import os
import subprocess
import time
from pathlib import Path
from collections import Counter

import pyrosetta
from dask.distributed import Client
from pyrosetta import *
from pyrosetta.distributed.cluster import PyRosettaCluster


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


def parse_mem_string(mem_string, multiply_by):
    # remove all letters from the string, save also the letters
    mem_int = int(''.join(filter(str.isdigit, mem_string)))
    print(mem_int)
    mem_unit = ''.join(filter(str.isalpha, mem_string))

    return str(mem_int * multiply_by) + mem_unit


def protocol(pose_path, min_rec_bb, index, output_path):
    """Run FlexPepDock protocol and write to worker-specific silent files."""
    from dask.distributed import get_worker
    import os
    import pyrosetta
    import pyrosetta.distributed.io as io
    import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job

    # Initialize worker-specific silent file
    worker = get_worker()
    worker_id = worker.id
    silent_file_name = os.path.join(output_path, f"decoys_worker_{worker_id}.silent")

    pose = io.pose_from_file(pose_path)
    basename = os.path.splitext(os.path.basename(pose_path))[0]
    numbered_basename = f"{basename}_{index:04d}"

    print(f'Running FlexPepDock on {numbered_basename} (Worker {worker_id})')

    # FlexPepDock XML definition
    xml = f"""
        <ROSETTASCRIPTS>
            <SCOREFXNS>
            </SCOREFXNS>
            <MOVERS>
                <FlexPepDock name="fpd" pep_refine="1" extra_scoring="1" lowres_preoptimize="1" min_receptor_bb='{int(min_rec_bb)}'/>
            </MOVERS>

            <PROTOCOLS>
                <Add mover_name="fpd" />
            </PROTOCOLS>
        </ROSETTASCRIPTS>
    """

    # Run RosettaScripts
    packed_pose = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)(pose.pose.clone())
    out_pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
    out_pose.pdb_info().name(numbered_basename)

    # Add extra scores to Pose object
    extra_scores = dict(get_string_real_pairs_from_current_job())
    for k, v in extra_scores.items():
        if not k.startswith('best'):
            out_pose.scores[k] = v

    # Save to the worker-specific silent file
    io.to_silent(out_pose, silent_file_name)
    print(f'Finished running FlexPepDock on {numbered_basename} (Worker {worker_id})')
    sys.stdout.flush()


def run_fpd_cluster(list_of_inputs, init_opts, min_rec_bb=False, output_path='.', nstruct=1):
    from dask.distributed import Client
    from dask_jobqueue import SLURMCluster
    import socket
    import os
    import sys

    pyrosetta.distributed.init(init_opts)

    # Setup the cluster
    tasks_per_job = os.environ['FPD_TASKS_PER_JOB']                                     # Number of tasks per node (one per core)
    memory_per_job = parse_mem_string(os.environ['FPD_MEM_PER_TASK'], tasks_per_job)    # Memory per task

    # Dashboard will be available at ip of localhost, port: 8787
    print(f"Dashboard URL: http://{socket.gethostname()}:8787")
    print(f"Submitting job array with {tasks_per_job} tasks per job, {memory_per_job} memory per job")
    sys.stdout.flush()
    sys.stderr.flush()

    cluster = SLURMCluster(
        cores=int(tasks_per_job),
        processes=int(tasks_per_job),
        memory=memory_per_job,
        walltime=os.environ['FPD_TIME'],
        interface="eth0",
        local_directory="/tmp",  # Temporary directory for workers
        name="fpd-${JOB_ID}",
        job_extra = [f"--array=1-{os.environ['FPD_NUM_JOBS']}"],
        env_extra=[
            "JOB_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}",
        ]
    )
    cluster.scale(jobs=1)
    client = Client(cluster)

    # Run tasks
    task_args = [(pose_path, min_rec_bb, i, output_path)
                 for pose_path in list_of_inputs for i in range(1, nstruct + 1)]

    futures = client.map(lambda args: protocol(*args), task_args)
    client.gather(futures)
    print("All tasks completed.")

    # Combine silent files, extract scores and save it to a separate file
    print('Combining silent files and extracting scores...')
    combined_silent_file = os.path.join(output_path, "decoys.silent")
    score_file_path = os.path.join(output_path, "score.sc")

    with open(combined_silent_file, "w") as combined_file, open(score_file_path, "w") as score_file:
        for silent_file in [os.path.join(output_path, f) for f in os.listdir(output_path) if f.startswith("decoys_worker_")]:
            with open(silent_file, "r") as infile:
                for line in infile:
                    if line.startswith("SCORE"):
                        score_file.write(line)
            os.remove(silent_file)  # Clean up individual worker silent files

    print(f"Combined silent file written to {combined_silent_file}")
    print(f"Scores written to {score_file_path}")


# Decide on nstruct based on upper limit or user input
def set_nstruct(n_templates, max_nstruct_per_decoy, max_nstruct_total=20000, args=None):
    if args and args.nstruct:
        nstruct = args.nstruct
    else:
        if max_nstruct_per_decoy * n_templates <= max_nstruct_total:
            nstruct = 10
        elif n_templates > max_nstruct_total:
            nstruct = 1
        else:
            nstruct = max_nstruct_total // n_templates # rounds down to closest integer

    print(f"Refinement will generate {nstruct} decoys per input structure for {n_templates} inputs, total {nstruct * n_templates} decoys")

    return nstruct


# Clean up existing output files before running FlexPepDock
def clean_before_run():
    for file in ["score.sc", "decoys.silent"]:
        Path(file).unlink(missing_ok=True)


# process redundancy and/or benchmarking
def filter_inputs(input_files, redundancy, benchmark):
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
def create_init_opts(args):
    init_opts = "-ex1 -ex2aro -use_input_sc -overwrite"
    # Add arguments to the dict
    if args.unboundrot:
        init_opts += f" -unboundrot {args.unboundrot}"
    if args.native:
        init_opts += f" -native {args.native}"

    return init_opts


def main():
    parser = argparse.ArgumentParser(description="Controls the FlexPepDock refinement part of the protocol.")
    parser.add_argument("-t", "--nstruct", type=int, help="Number of FlexPepDock decoys to generate. Default: 1")
    parser.add_argument("-m", "--min_rec_bb", action="store_true", help="Use receptor backbone minimization? Doubles runtime. Default: False")
    parser.add_argument("-u", "--unboundrot", type=str, help="Add unbound rotamers for the receptor. Default: None")
    parser.add_argument("-a", "--native", type=str, help="Native structure for comparison. Default: None")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Enable benchmarking mode")
    parser.add_argument("-r", "--redundancy", action="store_true", help="Filter for template redundancy. Default: False")
    parser.add_argument("-c", "--cpu", type=int, help="Number of CPUs to use. Default: All available")
    parser.add_argument("-o", "--outdir", type=str, help="Output directory. Default: Current directory", default=os.getcwd())
    args = parser.parse_args()

    # Collect input files
    input_files = glob.glob("*0001.pdb")

    # Filter input files
    input_files = filter_inputs(input_files, args.redundancy, args.benchmark)

    # Set number of decoys to generate
    nstruct = set_nstruct(len(input_files), 10, 20000, args)

    # For init Rosetta
    init_opts = create_init_opts(args)

    # Clean up existing output files before running FlexPepDock
    clean_before_run()

    # Open dask client and link to dashboard
    from dask.config import set
    set({'distributed.worker.memory.target': 0.8, 'distributed.worker.memory.spill': 0.75, 'distributed.worker.memory.pause':0.85})

    # Run FlexPepDock
    print('Running FlexPepDock...')
    run_fpd_cluster(input_files, init_opts, args.min_rec_bb, args.outdir, int(nstruct))


if __name__ == "__main__":
    from dask.distributed import Client
    main()


