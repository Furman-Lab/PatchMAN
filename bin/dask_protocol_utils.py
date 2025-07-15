import os
import subprocess
from dask.distributed import Client, LocalCluster, as_completed
from datetime import datetime

cluster = None  # Global variable to hold the Dask cluster instance

def run_dask_tasks(cluster, task_func, task_args_list, task_name):
	"""Execute tasks using Dask cluster and monitor progress"""
	client = Client(cluster.scheduler_address)

	print(f"Starting {task_name} with {len(task_args_list)} tasks")
	print(f"Dask dashboard: {client.dashboard_link}")

	# Submit all tasks
	futures = []
	for args in task_args_list:
		future = client.submit(task_func, *args)
		futures.append(future)

	# Process results as they complete
	completed_count = 0
	failed_tasks = []

	for future in as_completed(futures):
		try:
			result = future.result()
			completed_count += 1

			# Create individual log file for this task
			task_id = result.get('task_id', 'unknown')
			# Clean task_id for filename (remove problematic characters)
			clean_task_id = str(task_id).replace('/', '_').replace('\\', '_').replace(':', '_')
			log_filename = f"{task_name}_{clean_task_id}.log"

			with open(log_filename, "w") as task_log:
				task_log.write(f"Task Name: {task_name}\n")
				task_log.write(f"Task ID: {task_id}\n")
				task_log.write(f"Return Code: {result['returncode']}\n")
				task_log.write(f"Timestamp: {datetime.now().isoformat()}\n")
				task_log.write("=" * 60 + "\n\n")

				if result.get('stdout'):
					task_log.write("STDOUT:\n")
					task_log.write("-" * 40 + "\n")
					task_log.write(result['stdout'])
					task_log.write("\n" + "-" * 40 + "\n\n")

				if result.get('stderr'):
					task_log.write("STDERR:\n")
					task_log.write("-" * 40 + "\n")
					task_log.write(result['stderr'])
					task_log.write("\n" + "-" * 40 + "\n\n")

			if result['returncode'] != 0:
				failed_tasks.append(result)
				print(f"Task {result['task_id']} failed with return code {result['returncode']} - Log: {log_filename}")
				if result['stderr']:
					print(f"Error: {result['stderr']}")
			else:
				print(
					f"Task {result['task_id']} completed successfully ({completed_count}/{len(task_args_list)}) - Log: {log_filename}")

		except Exception as e:
			print(f"Task failed with exception: {e}")
			failed_tasks.append({'task_id': 'unknown', 'error': str(e)})

			# Create log file for exception
			exception_log_filename = f"{task_name}_exception_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
			with open(exception_log_filename, "w") as exc_log:
				exc_log.write(f"Task Name: {task_name}\n")
				exc_log.write(f"Exception: {str(e)}\n")
				exc_log.write(f"Timestamp: {datetime.now().isoformat()}\n")

			print(f"Exception logged to: {exception_log_filename}")

	client.close()

	if failed_tasks:
		print(f"Warning: {len(failed_tasks)} tasks failed in {task_name}")
		return False
	else:
		print(f"All {task_name} tasks completed successfully")
		return True

def setup_cluster(n_searches, cpu_count):
	"""Get existing cluster or create new one if needed"""
	global cluster
	if cluster is None:
		cluster = LocalCluster(
			n_workers=min(cpu_count, n_searches),
			threads_per_worker=1,
			memory_limit='500MB',
			dashboard_address=':8787'
		)
		print(f"Created Dask cluster with {cluster.workers} workers")
	return cluster


def run_master_task(task_id, master_cutoff, master_args):
	"""Execute MASTER task for a single motif"""
	import glob

	# Get list of .pds files (equivalent to master_list=($(ls *pds)))
	master_list = glob.glob("*.pds")
	master_list.sort()  # Ensure consistent ordering

	if task_id >= len(master_list):
		return {
			'task_id': task_id,
			'returncode': 1,
			'stdout': '',
			'stderr': f'Task ID {task_id} exceeds number of .pds files ({len(master_list)})'
		}

	# Get the specific .pds file for this task
	master_job = master_list[task_id]
	print(task_id, master_job)

	print(f"Task {task_id}: Processing {master_job}")

	# Run MASTER command directly (equivalent to the original master command)
	cmd = ['master', '--query', master_job, '--bbRMSD', '--rmsdCut', str(master_cutoff),
		   '--targetList', master_args, '--topN', '100000', '--matchOut', f'{master_job}.matches']

	result = subprocess.run(cmd, capture_output=True, text=True)
	return {
		'task_id': task_id,
		'returncode': result.returncode,
		'stdout': result.stdout,
		'stderr': result.stderr,
		'master_job': master_job
	}


def extract_templates_task(task_id, ppkrec, focus_mask_args):
	"""Execute template extraction task for a single motif"""
	import glob

	# Get list of .matches files (equivalent to matches=($(ls *matches)))
	matches = glob.glob("*.matches")
	matches.sort()  # Ensure consistent ordering

	if task_id >= len(matches):
		return {
			'task_id': task_id,
			'returncode': 1,
			'stdout': '',
			'stderr': f'Task ID {task_id} exceeds number of .matches files ({len(matches)})'
		}

	# Get the specific .matches file for this task
	match_list = matches[task_id]

	# Read motif_list file to get the corresponding motif
	with open("motif_list", "r") as f:
		motifs = [line.strip() for line in f.readlines()]

	if task_id >= len(motifs):
		return {
			'task_id': task_id,
			'returncode': 1,
			'stdout': '',
			'stderr': f'Task ID {task_id} exceeds number of motifs ({len(motifs)})'
		}

	motif = motifs[task_id]

	print(f"Task {task_id}: Processing {match_list} with motif {motif}")

	# Set environment variable
	env = os.environ.copy()
	env['PYTHONNOUSERSITE'] = '1'

	# Run the extract_peps_for_motif.py script
	cmd = ['python3', f"{os.environ['BIN_DIR']}/extract_peps_for_motif.py",
		   '-m', match_list, '-p', 'pepfile', '-r', ppkrec, '--patch', motif] + focus_mask_args

	result = subprocess.run(cmd, env=env, capture_output=True, text=True)

	return {
		'task_id': task_id,
		'returncode': result.returncode,
		'stdout': result.stdout,
		'stderr': result.stderr,
		'match_file': match_list,
		'motif': motif
	}
