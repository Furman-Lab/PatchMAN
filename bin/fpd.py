import argparse
import glob
import os
import socket
import sys
from collections import Counter
from pathlib import Path

import pyrosetta
from dask_jobqueue import SLURMCluster
from pyrosetta import *
from dask.distributed import Client, WorkerPlugin


class PyRosettaFPDPlugin(WorkerPlugin):
    def __init__(self, init_options, min_rec_bb=False, native_path=None):
        self.init_options = init_options
        self.min_rec_bb = min_rec_bb
        self.native_path = native_path  # Path to native PDB if needed
    
    def setup(self, worker):
        """Initialized FPD mover only once for each worker"""
        import pyrosetta
        from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
        self.worker = worker
        
        # Initialize PyRosetta once
        print(f"Initializing PyRosetta on worker {worker.id}")
        pyrosetta.init(self.init_options)
        
        # Create the FPD mover once
        print(f"Creating FlexPepDock mover on worker {worker.id}")

        # FlexPepDock XML definition
        xml = f"""
            <ROSETTASCRIPTS>
                <SCOREFXNS>
                    <ScoreFunction name="sfxn" weights="ref2015"/>
                </SCOREFXNS>
                <MOVERS>
                    <FlexPepDock name="fpd" pep_refine="1" extra_scoring="1" lowres_preoptimize="1" min_receptor_bb='{int(self.min_rec_bb)}'/>
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

        # Set the native pose if provided
        if self.native_path:
            native_pose = io.pose_from_file(self.native_path)
            fpd_protocol.set_native_pose(native_pose)
            print(f"Native pose set on worker {worker.id} from {self.native_path}")
        else:
            print("No native pose provided, using the input as native.")

        # Store the mover as a worker attribute
        worker.fpd_mover = fpd_protocol

        print(f"PyRosetta and FlexPepDock initialized on worker {worker.id}")
        sys.stdout.flush()


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
    """ We need to calculate the required memory per task based on the number of tasks per job """
    
    # remove all letters from the string, get also the letters separately
    mem_int = int(''.join(filter(str.isdigit, mem_string)))
    mem_unit = ''.join(filter(str.isalpha, mem_string))
    
    return str(mem_int * multiply_by) + mem_unit


def get_cluster_args(list_of_inputs, nstruct):
    """Get the argument values from the command line"""
    # Set up the cluster
    if 'FPD_NUM_JOBS' not in os.environ:
        os.environ['FPD_TASKS_PER_JOB'] = str(os.cpu_count() // 2)
    tasks_per_job = int(os.environ['FPD_TASKS_PER_JOB'])  # Number of tasks per node (one per core)
    
    if 'FPD_MEM_PER_TASK' not in os.environ:
        memory_per_job = f'{tasks_per_job * 2}G'
    else:
        memory_per_job = parse_mem_string(os.environ['FPD_MEM_PER_TASK'], tasks_per_job)  # Memory per task
    
    if 'FPD_TIME' not in os.environ:
        os.environ['FPD_TIME'] = "24:00:00"
    
    if 'FPD_NUM_JOBS' not in os.environ:
        os.environ['FPD_NUM_JOBS'] = str(len(list_of_inputs) * nstruct)
        
    return tasks_per_job, memory_per_job


def start_cluster(tasks_per_job, memory_per_job):
    """
    Start a Dask cluster on the SLURM scheduler and return a client object.
    """
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
        job_extra_directives = [f"--array=1-{os.environ['FPD_NUM_JOBS']}"],
        job_script_prologue=[
            "JOB_ID=${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}",
        ]
    )
    cluster.scale(jobs=1)
    client = Client(cluster)
    
    return client


def write_score_file_from_silent(output_path):
    """
    Combine silent files, extract scores and save it to a separate file
    """
    
    print('Combining silent files and extracting scores...')
    combined_silent_file = os.path.join(output_path, "decoys.silent")
    score_file_path = os.path.join(output_path, "score.sc")

    with open(combined_silent_file, "w") as combined_file, open(score_file_path, "w") as score_file:
        for silent_file in [os.path.join(output_path, f) for f in os.listdir(output_path) if f.startswith("decoys_worker_")]:
            with open(silent_file, "r") as infile:
                for line in infile:
                    if line.startswith("SCORE"):
                        score_file.write(line)
                    combined_file.write(line)
            os.remove(silent_file)  # Clean up individual worker silent files

    print(f"Combined silent file written to {combined_silent_file}")
    print(f"Scores written to {score_file_path}")
    
    
def create_task_list(list_of_inputs, nstruct=1, output_path='.', finished_decoys=None):
    """
    Create a list of tasks to run on the cluster. If finished_decoys is provided, skip those decoys.
    """
    # Create the task list
    if finished_decoys is None:
        task_args = [(pose_path, i, output_path)
                 for pose_path in list_of_inputs for i in range(1, nstruct + 1)]
        print(f"Running FlexPepDock on {len(list_of_inputs) * nstruct} decoys")
    else:
        task_args = []
        for pose_path in list_of_inputs:
            basename = os.path.splitext(os.path.basename(pose_path))[0]
            for i in range(1, nstruct + 1):
                if f"{basename}_{i:04d}" not in finished_decoys:
                    task_args.append((pose_path, i, output_path))
        print(f"Resuming with {len(task_args)} decoys to run out of {len(list_of_inputs) * nstruct} decoys")
        
    return task_args


def protocol(pose_path, index, output_path):
    """
    Run FlexPepDock protocol and write to worker-specific silent files.
    """
    from dask.distributed import get_worker
    import os
    import pyrosetta
    import pyrosetta.distributed.io as io
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job

    # Initialize worker-specific silent file
    worker = get_worker()
    worker_id = worker.id
    silent_file_name = os.path.join(output_path, f"decoys_worker_{worker_id}.silent")

    pose = io.pose_from_file(pose_path)
    basename = os.path.splitext(os.path.basename(pose_path))[0]
    numbered_basename = f"{basename}_{index:04d}"

    print(f'Running FlexPepDock on {numbered_basename} (Worker {worker_id})')

    out_pose = pyrosetta.distributed.packed_pose.to_pose(pose)
    fpd_protocol = worker.fpd_mover
    fpd_protocol.apply(out_pose)
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


def run_fpd_cluster(list_of_inputs, init_opts, min_rec_bb=False, output_path='.', nstruct=1, finished_decoys=None, native=None):
    """
    Initializes the cluster and its workers, generate the list of tasks to run.
    Then, combines the silent files and creates a score file.
    """
    pyrosetta.distributed.init(init_opts)
    
    # Setup the cluster
    tasks_per_job, memory_per_job = get_cluster_args(list_of_inputs, nstruct)
    
    # Create the task list
    task_args = create_task_list(list_of_inputs, nstruct, output_path, finished_decoys)
    
    # Run the cluster and the jobs
    client = start_cluster(tasks_per_job, memory_per_job)
    client.register_plugin(PyRosettaFPDPlugin(init_opts, min_rec_bb, native))
    futures = client.map(lambda args: protocol(*args), task_args)
    client.gather(futures)
    
    print("All tasks completed. Combining silent files and extracting scores...")
    write_score_file_from_silent(output_path)

# Decide on nstruct based on upper limit or user input
def set_nstruct(n_templates, max_nstruct_per_decoy, max_nstruct_total=20000, args=None):
    """
    User input always overwrites calcilation.
    Otherwise, we do not want to refine more then 20000 decoys, but we also dont want to refine more than 10 decoys per input structure.
    """
    if args and args.nstruct:
        nstruct = args.nstruct
    else:
        if max_nstruct_per_decoy * n_templates <= max_nstruct_total:
            nstruct = max_nstruct_per_decoy
        elif n_templates > max_nstruct_total:
            nstruct = 1
        else:
            nstruct = max_nstruct_total // n_templates # rounds down to the closest integer

    print(f"Refinement will generate {nstruct} decoys per input structure for {n_templates} inputs, total {nstruct * n_templates} decoys")

    return nstruct


# Clean up existing output files before running FlexPepDock
def clean_before_run(force_rerun):
    for file in ["score.sc", "decoys.silent"]:
        Path(file).unlink(missing_ok=True)
        
    if force_rerun:
        print("Deleting half-finished jobs...")
        for file in glob.glob("decoys_worker_*.silent"):
            Path(file).unlink()


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

    return init_opts


def list_decoys_for_restarting():
    # We need the list of decoys that were already run to restart the process
    finished_decoys = []
    for file in glob.glob("*.silent"):
        with open(file, "r") as infile:
            for line in infile:
                if line.startswith("SCORE") and "description" not in line:
                    finished_decoys.append(line.split()[-1])
    
    return set(finished_decoys)


def main():
    parser = argparse.ArgumentParser(description="Controls the FlexPepDock refinement part of the protocol.")
    parser.add_argument("-t", "--nstruct", type=int, help="Number of FlexPepDock decoys to generate. Default: 1")
    parser.add_argument("-m", "--min_rec_bb", action="store_true", help="Use receptor backbone minimization? Doubles runtime. Default: False")
    parser.add_argument("-u", "--unboundrot", type=str, help="Add unbound rotamers for the receptor. Default: None")
    parser.add_argument("-a", "--native", type=str, help="Native structure for comparison. Default: None")
    parser.add_argument("-b", "--benchmark", action="store_true", help="Enable benchmarking mode")
    parser.add_argument("-f", "--force_rerun", action="store_true", help="Force rerunning jobs that stopped in the mid. Default: False")
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
    clean_before_run(args.force_rerun)
    
    # If we did not want to force rerunning and there are some decoys already finished, we need to list them
    finished_decoys = list_decoys_for_restarting() if not args.force_rerun else None
    
    # If args.native specified, check if it exists. Exit if not
    if args.native and not Path(args.native).exists():
        raise FileNotFoundError(f"Native structure {args.native} not found. Exiting...")
    
    # Run FlexPepDock
    print('Running FlexPepDock...')
    run_fpd_cluster(input_files, init_opts, args.min_rec_bb, args.outdir, int(nstruct), finished_decoys, args.native)


if __name__ == "__main__":
    from dask.distributed import Client
    main()


