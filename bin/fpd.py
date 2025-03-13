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
        import pyrosetta
        
        # Step 1: Initialize PyRosetta once
        print(f"Initializing PyRosetta on worker {worker.id}")
        pyrosetta.init(self.init_options)
        
        # Step 2: Load native pose if path is provided
        native_pose = None
        if self.native_path:
            print(f"Loading native pose from {self.native_path}")
            native_pose = pyrosetta.pose_from_pdb(self.native_path)
        
        # Step 3: Create the FPD mover once
        print(f"Creating FlexPepDock mover on worker {worker.id}")
        
        from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
        
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
        if native_pose:
            fpd_protocol.set_native_pose(native_pose)
            print(f"Native pose set on worker {worker.id}")
        else:
            print("No native pose provided, using the input as native.")
        
        # Store the mover as a worker attribute for reuse
        worker.fpd_mover = fpd_protocol


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
    # remove all letters from the string, get also the letters separately
    mem_int = int(''.join(filter(str.isdigit, mem_string)))
    mem_unit = ''.join(filter(str.isalpha, mem_string))
    
    return str(mem_int * multiply_by) + mem_unit

def create_fpd_mover(native_pose=None, min_rec_bb=False):
    from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
    # FlexPepDock XML definition
    xml = f"""
        <ROSETTASCRIPTS>
            <SCOREFXNS>
                <ScoreFunction name="sfxn" weights="ref2015"/>
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
    parser = RosettaScriptsParser()
    default_options = pyrosetta.rosetta.basic.options.process()
    tag = parser.create_tag_from_xml_string(xml, default_options)
    fpd_protocol = parser.parse_protocol_tag(tag, default_options).get_mover(1)
    fpd_protocol.set_default() # this will take care of unboundrot
    
    # Set the native pose if provided
    if native_pose:
        fpd_protocol.set_native_pose(native_pose)
        print(fpd_protocol.get_native_pose().pdb_info())
    else:
        print("No native pose provided, using the input as native.")
        
    return fpd_protocol


def protocol(pose_path, min_rec_bb, index, output_path, native_pose=None):
    """Run FlexPepDock protocol and write to worker-specific silent files."""
    from dask.distributed import get_worker
    import os
    import pyrosetta
    import pyrosetta.distributed.io as io
    from pyrosetta.rosetta.protocols.rosetta_scripts import RosettaScriptsParser
    from pyrosetta.rosetta.protocols.jd2 import get_string_real_pairs_from_current_job

    # Initialize worker-specific silent file
    worker = get_worker()
    worker_id = worker.id
    silent_file_name = os.path.join(output_path, f"decoys_worker_{worker_id}.silent")

    pose = io.pose_from_file(pose_path)
    basename = os.path.splitext(os.path.basename(pose_path))[0]
    numbered_basename = f"{basename}_{index:04d}"
    
    default_options = pyrosetta.rosetta.basic.options.process()
    print(f'Parsed command line options: {default_options.get_argv()}')

    print(f'Running FlexPepDock on {numbered_basename} (Worker {worker_id})')

    # FlexPepDock XML definition
    xml = f"""
        <ROSETTASCRIPTS>
            <SCOREFXNS>
                <ScoreFunction name="sfxn" weights="ref2015"/>
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
    parser = RosettaScriptsParser()
    default_options = pyrosetta.rosetta.basic.options.process()
    tag = parser.create_tag_from_xml_string(xml, default_options)
    fpd_protocol = parser.parse_protocol_tag(tag, default_options).get_mover(1)
    fpd_protocol.set_default() # this will take care of unboundrot
    
    # Set the native pose if provided
    if native_pose:
        fpd_protocol.set_native_pose(native_pose)
        print(fpd_protocol.get_native_pose().pdb_info())
    else:
        print("No native pose provided, using the input as native.")
        
    # Run the protocol on a copy of the input pose
    out_pose = pose.pose.clone()
    fpd_protocol.apply(out_pose)
    
    #packed_pose = rosetta_scripts.SingleoutputRosettaScriptsTask(xml)(pose.pose.clone())
    #out_pose = pyrosetta.distributed.packed_pose.to_pose(packed_pose)
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
    
    
def setup_cluster(list_of_inputs, nstruct):
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


def start_cluster(tasks_per_job, memory_per_job, init_opts):
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
    client.register_plugin(PyRosettaInitPlugin(init_opts))
    
    return client


def write_score_file_from_silent(output_path):
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
                    combined_file.write(line)
            os.remove(silent_file)  # Clean up individual worker silent files

    print(f"Combined silent file written to {combined_silent_file}")
    print(f"Scores written to {score_file_path}")
    
    
def create_task_list(list_of_inputs, fpd_mover, nstruct=1, output_path='.', finished_decoys=None):
    # Create the task list
    if finished_decoys is None:
        task_args = [(pose_path, fpd_mover, i, output_path)
                 for pose_path in list_of_inputs for i in range(1, nstruct + 1)]
        print(f"Running FlexPepDock on {len(list_of_inputs) * nstruct} decoys")
    else:
        task_args = []
        for pose_path in list_of_inputs:
            basename = os.path.splitext(os.path.basename(pose_path))[0]
            for i in range(1, nstruct + 1):
                if f"{basename}_{i:04d}" not in finished_decoys:
                    task_args.append((pose_path, fpd_mover, i, output_path))
        print(f"Resuming with {len(task_args)} decoys to run out of {len(list_of_inputs) * nstruct} decoys")
        
    return task_args


def run_fpd_cluster(list_of_inputs, init_opts, min_rec_bb=False, output_path='.', nstruct=1, finished_decoys=None, native=None):
    """
    Run FlexPepDock protocol on a list of input structures using a Dask cluster.
    :param list_of_inputs: List of input structures, as paths
    :param init_opts: Options to initialize PyRosetta, a string
    :param min_rec_bb: Use receptor backbone minimization
    :param output_path: Output directory
    :param nstruct: Number of decoys to generate per input structure
    :param finished_decoys: List of decoys that were already run, so that we not repeat them
    :param native: Native structure for comparison
    """
    
    print(f'starting pyrosetta with {init_opts}' )
    pyrosetta.init(init_opts)
    
    # Set up the cluster
    tasks_per_job, memory_per_job = setup_cluster(list_of_inputs, nstruct)
    
    # Create the task list
    native_pose = io.pose_from_file(native) if native else None
    fpd_mover = create_fpd_mover(native_pose, min_rec_bb)
    task_args = create_task_list(list_of_inputs, fpd_mover, nstruct, output_path, finished_decoys)
    
    # Start the cluster
    client = start_cluster(tasks_per_job, memory_per_job, init_opts)
    futures = client.map(lambda args: protocol(*args), task_args)
    client.gather(futures)
    print("All tasks completed.")

    # Extract the scores from the silent file to a separate file
    write_score_file_from_silent(output_path)


# Decide on nstruct based on upper limit or user input
def set_nstruct(n_templates, max_nstruct_per_decoy, max_nstruct_total=20000, args=None):
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
    
    # Run FlexPepDock
    print('Running FlexPepDock...')
    run_fpd_cluster(input_files, init_opts, args.min_rec_bb, args.outdir, int(nstruct), finished_decoys, args.native)


if __name__ == "__main__":
    from dask.distributed import Client
    main()


