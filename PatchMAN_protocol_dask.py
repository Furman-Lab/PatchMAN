#!/usr/bin/env python3
from bin.protocol_utils import *
from bin.dask_protocol_utils import *
from dask.distributed import Client, LocalCluster, as_completed
from dask import delayed
import subprocess
import os

def main():
	# Gather all the arguments needed for the protocol
	args = parse_protocol_args()
	config = load_config('config.ini')
	os.chdir(args.work_dir)
	cluster = None

	# Prepare inputs
	pep_sequence = validate_peptide_sequence(args.peptide_sequence)
	args.focus_mask_args = prepare_mask_and_focus(args)

	# Prepare receptor and derived filenames 
	receptor, receptor_base, ppkrec = prepare_and_set_filenames(args.receptor)
	args.fpd_args.append(f"-u {receptor}")
	args.fpd_args.append(f"--use_local")

	# Copy also native file into directory
	if args.native_pdb:
		args.native_pdb = prepare_pdb(args.native_pdb)

	# Prepare for benchmark mode if provided
	args = prepare_benchmark_mode(args, receptor_base)

	# Step 1: Split to motifs
	pdb_list_file = "motif_list"
	if 1 in args.steps_range:
		print("Step 1: Splitting to motifs...")
		results_split_to_motifs = subprocess.run(['python3', f"{os.environ['BIN_DIR']}/split_to_motifs.py", "-i", receptor,
						*args.focus_mask_args.split()], check=True)
		
		check_subprocess_result(results_split_to_motifs, 1)
	
	# in hotspot mode, Step 1 also created a new .pdb called focus_from_hotspots.pdb - this needs to be changed for the extraction step
	if args.hotspot_mode:
		args.focus_mask_args = "-o -s focus_from_hotspots.pdb"

	# Calculate n_searches using motif_list
	n_searches = sum(1 for line in open(pdb_list_file))
	print(f"Number of searches: {n_searches}")

	# Step 2: Prepare input (run directly)
	if 2 in args.steps_range:
		print("Step 2: Preparing input...")

		clean_pdbs = glob.glob("*.clean.pdb")
		if clean_pdbs:
			result_prepack = subprocess.run(['python3', f"{os.environ['BIN_DIR']}/prepack_receptor.py", clean_pdbs[0]],
									capture_output=True, text=True)
			
			check_subprocess_result(result_prepack, 2)
		else:
			raise FileNotFoundError(
				"No .clean.pdb files found in the working directory. Probably there was a problem with generating them.")

	# Step 3: Run MASTER
	if 3 in args.steps_range:
		print("Step 3: Running MASTER search...")
		cluster = setup_cluster(n_searches, args.cpu)

		# Prepare task arguments for all motifs
		master_task_args = [(i, args.master_cutoff, args.master_args)
							for i in range(n_searches)]

		master_success = run_dask_tasks(
			cluster, run_master_task, master_task_args, "master_search"
		)
		check_dask_tasks(master_success, 3)

	# Step 4: Extract templates
	if 4 in args.steps_range:
		print("Step 4: Extracting templates...")
		cluster = setup_cluster(n_searches, args.cpu)

		# Create pepfile with peptide sequence
		pepfile = f"pepfile"
		with open(pepfile, "w") as f:
			f.write(args.peptide_sequence)

		# Prepare task arguments for template extraction
		extract_task_args = [(i, ppkrec, args.focus_mask_args.split())
							 for i in range(n_searches)]

		extract_success = run_dask_tasks(
			cluster, extract_templates_task, extract_task_args, "template_extraction"
		)
		
		check_dask_tasks(extract_success, 4)

	# Clean up Dask cluster before Step 5 to free up resources
	if cluster:
		cluster.close()
		print("Dask cluster closed")

	# Step 5: FlexPepDock (FPD) - Run directly, no Slurm submission
	if 5 in args.steps_range:
		print("Step 5: Running FlexPepDock...")
		print(f"FPD called with arguments: {' '.join(args.fpd_args)}")

		# Run FlexPepDock Python script directly (using PROTOCOL_ROOT path like original)
		fpd_cmd = ['python3', f"{os.environ['PROTOCOL_ROOT']}/bin/fpd.py"] + args.fpd_args

		# Filter out warnings like the original script
		for line in fpd_process.stdout:
			if 'WARNING' not in line.upper():
				print(line.rstrip())
		
		# Then check the result at the end
		fpd_success = fpd_process.wait()
		if fpd_success != 0:
			print(f"Step {step_number} failed with return code {return_code}")
			raise subprocess.CalledProcessError(return_code, fpd_process.args)
		else:
			print(f"Step {step_number} completed successfully")

	# Step 6: Clustering & Finalizing - Run directly, no Slurm submission
	if 6 in args.steps_range:
		print("Step 6: Clustering and finalizing...")

		# Create clustering directory (equivalent to mkdir -p clustering/)
		os.makedirs("clustering", exist_ok=True)

		# Run clustering Python script and redirect output to clog file
		cluster_cmd = ['python3', f"{os.environ['PROTOCOL_ROOT']}/bin/cluster.py"]

		with open("clustering/clog", "w") as clog_file:
			cluster_result = subprocess.run(cluster_cmd, stdout=clog_file, stderr=subprocess.STDOUT, text=True)
			
			check_subprocess_result(result_prepack, 2)

	# Print summary
	print(f'\nCompleted processing of receptor {receptor_base} with peptide: {args.peptide_sequence}')
	print("Execution summary:")
	print(f"  Step 1 (Split motifs): {'✓ Completed' if 1 in args.steps_range else '- Skipped'}")
	print(f"  Step 2 (Prepare input): {'✓ Completed' if 2 in args.steps_range else '- Skipped'}")
	print(f"  Step 3 (MASTER): {'✓ Completed with Dask' if 3 in args.steps_range else '- Skipped'}")
	print(f"  Step 4 (Extract templates): {'✓ Completed with Dask' if 4 in args.steps_range else '- Skipped'}")
	print(f"  Step 5 (FlexPepDock): {'✓ Completed' if 5 in args.steps_range else '- Skipped'}")
	print(f"  Step 6 (Clustering): {'✓ Completed' if 6 in args.steps_range else '- Skipped'}")

	print("\nAll steps completed successfully!")


if __name__ == '__main__':
	main()
