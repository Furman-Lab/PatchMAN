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


def protocol(pose_path, min_rec_bb, index):
    import pyrosetta
    import pyrosetta.distributed.io as io
    import pyrosetta.distributed.tasks.rosetta_scripts as rosetta_scripts
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job

    pose = io.pose_from_file(pose_path)
    basename = os.path.splitext(os.path.basename(pose_path))[0]
    numbered_basename = f"{basename}_{index:04d}"

    print(f'Running FlexPepDock on {numbered_basename}')

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

    # Name the decoys according to the usual Rosetta naming scheme. Need to convert to regular Pose
    packed_pose = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)(pose.pose.clone())
    out_pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
    out_pose.pdb_info().name(numbered_basename)

    # Add extra scores to Pose object
    extra_scores = dict(get_string_real_pairs_from_current_job())
    for k, v in extra_scores.items():
        if not k.startswith('best'):
            out_pose.scores[k] = v

    # Save into the silent file
    io.to_silent(out_pose, "decoys.silent")
    print(f'Finished running FlexPepDock on {numbered_basename}')
    sys.stdout.flush()
    sys.stderr.flush()


def run_fpd_cluster2(list_of_inputs, init_opts, min_rec_bb=False, output_path='.', nstruct=1):
    from dask.distributed import Client, LocalCluster
    import socket
    pyrosetta.distributed.init(init_opts)

    # Set up a LocalCluster using resources provided by SLURM
    # If the cluster starts to have 'Too many files open', increase the number of tasks per worker
    task_per_worker = 2
    cores_per_worker = int(os.getenv("SLURM_CPUS_PER_TASK", 1)) * task_per_worker
    total_tasks = int(int(os.getenv("SLURM_NTASKS", 100))/task_per_worker)  
    memory_limit = f'{str(int(os.getenv("SLURM_MEM_PER_CPU", "1600")) * task_per_worker)}MB' 
    host_ip = socket.gethostbyname(socket.gethostname())

    print(f"Starting LocalCluster with {total_tasks} workers, {cores_per_worker} cores per worker, {memory_limit} memory per worker on {host_ip}")
    sys.stdout.flush()
    cluster = LocalCluster(n_workers=total_tasks,
                           threads_per_worker=cores_per_worker,
                           memory_limit=memory_limit,
                           local_directory='/tmp',
                           processes=True,
                           dashboard_address=f'{host_ip}:8787',
                           host=host_ip
                           )
    client = Client(cluster)
    print(f"Dashboard URL: {client.dashboard_link}")
    sys.stdout.flush()
    sys.stderr.flush()

    # Run jobs
    futures = [client.submit(protocol, pose_path, min_rec_bb, i)
           for pose_path in list_of_inputs for i in range(1, nstruct + 1)]

    print(f'FUTURES: {len(futures)}')
    results = client.gather(futures)

    # Make a score file: open the silent file and extract lines starting with 'SCORE'
    print('FPD run finished, collecting scores from silent file')
    with open(os.path.join(output_path, "decoys.silent"), "r") as infile:
        score_lines = [line for line in infile if line.startswith("SCORE")]

    # Write the extracted lines to a new file
    print('Writing scores to file')
    with open(os.path.join(output_path, "score.sc"), "w") as outfile:
        outfile.writelines(score_lines)


def run_fpd_cluster(list_of_inputs, init_opts, min_rec_bb=False, output_path='.', nstruct=1):
    from joblib import Parallel, delayed
    import os
    import socket
    import pyrosetta.distributed.io as io
    from pyrosetta import init

    # Initialize PyRosetta
    init(init_opts)

    # Split workload across workers
    total_tasks = len(list_of_inputs) * nstruct
    n_jobs = int(os.getenv("SLURM_NTASKS", 1))
    print(f"Starting parallel execution with {total_tasks} tasks...")
    sys.stdout.flush()

    # Generate task list
    tasks = [
        delayed(protocol)(pose_path, min_rec_bb, i)
        for pose_path in list_of_inputs
        for i in range(1, nstruct + 1)
    ]

    # Execute tasks in parallel
    results = Parallel(n_jobs=n_jobs, backend="loky")(tasks)

    # Process results
    print(f"Collected {len(results)} results. Writing scores...")
    sys.stdout.flush()

    # Extract scores from the silent file
    score_file_path = os.path.join(output_path, "score.sc")
    with open(os.path.join(output_path, "decoys.silent"), "r") as infile:
        score_lines = [line for line in infile if line.startswith("SCORE")]

    with open(score_file_path, "w") as outfile:
        outfile.writelines(score_lines)

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


# Calculate the number of processors to use
def set_procs(args, nstruct, input_files):
    cpu_count = os.cpu_count() - 1
    if args.cpu:
        print('Using user-defined CPU count: ', args.cpu)
        procs = min(args.cpu - 1, cpu_count)
    elif "SLURM_NTASKS" not in os.environ:
        print('Using all available CPUs, SLURM_NTASKS is not available: ', cpu_count)
        procs = min(nstruct * len(input_files), cpu_count)
    else:
        print('Using SLURM_NTASKS: ', os.environ["SLURM_NTASKS"])
        procs = int(os.environ["SLURM_NTASKS"]) - 1

    print(f"Running FlexPepDock refinement on {procs} threads")

    return procs


# Clean up existing output files before running FlexPepDock
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
    input_files = filter(input_files, args.redundancy, args.benchmark)

    # Set number of decoys to generate
    nstruct = set_nstruct(len(input_files), 10, 20000, args)

    # Set number of processors to use
    procs = set_procs(args, nstruct, input_files)

    # For init Rosetta
    init_opts = create_init_opts(args)

    # Clean up existing output files before running FlexPepDock
    clean_before_run()

    # Open dask client and link to dashboard
    procs = int(procs / 2) # for debugging
    from dask.config import set
    set({'distributed.worker.memory.target': 0.8, 'distributed.worker.memory.spill': 0.75, 'distributed.worker.memory.pause':0.85})

    # Run FlexPepDock
    run_fpd_cluster2(input_files, init_opts, args.min_rec_bb, args.outdir, int(nstruct))


if __name__ == "__main__":
    from dask.distributed import Client
    main()


