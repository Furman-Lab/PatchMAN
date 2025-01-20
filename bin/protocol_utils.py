import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
import configparser
import socket
import glob

def load_config(config_file_path="config.ini"):
	'''
	Load the configuration file and set up environment variables
	:param config_file_path: Path to the configuration file
	:return: Dictionary containing the configuration settings
	'''
	config = configparser.ConfigParser()
	config.read(config_file_path)
	
	config_dict = dict()
	
	# Determine the hostname section to use
	hostname = socket.gethostname()
	specific_config = ''
	for section in config.sections():
		if section.startswith('HOSTNAME:') and section[len('HOSTNAME:'):] in hostname:
			specific_config = section
			break

	# Handle specific configurations based on hostname
	config_dict['ADD_SBATCH'] = config.get(specific_config, 'ADD_SBATCH', fallback='')
	os.environ['FPD_NUM_JOBS'] = str(config.get(specific_config, 'FPD_NUM_JOBS', fallback=50))
	os.environ['FPD_TASKS_PER_JOB'] = str(config.get(specific_config, 'FPD_TASKS_PER_JOB', fallback=3))
	os.environ['FPD_MEM_PER_TASK'] = str(config.get(specific_config, 'FPD_MEM_PER_TASK', fallback='2G'))
	os.environ['FPD_TIME'] = str(config.get(specific_config, 'FPD_TIME', fallback='120:00:00'))

	# Initialize PROTOCOL_ROOT from the configuration or use the script's location as a default
	protocol_root = Path(config.get('DEFAULT', 'PROTOCOL_ROOT', fallback=str(Path(__file__).resolve().parent.parent)))
	
	# Dynamically derive DB_PATH from PROTOCOL_ROOT if not specified
	db_path = config.get('DEFAULT', 'DB_PATH', fallback=protocol_root / 'databases/master_clean/')
	
	# Setting up environment variables
	os.environ['PROTOCOL_ROOT'] = str(protocol_root)
	os.environ['DB_PATH'] = str(db_path)
	os.environ['PYTHON'] = 'python3'
	os.environ['BIN_DIR'] = f'{protocol_root}/bin'
	
	# Example: Print environment variables to verify
	print("Environment Variables Set:")
	print(f"PROTOCOL_ROOT={os.environ['PROTOCOL_ROOT']}")
	print(f"DB_PATH={os.environ['DB_PATH']}")
	print(f"ADD_SBATCH={config_dict['ADD_SBATCH']}")

	# check if sbatch command is available
	try:
		subprocess.run(["sbatch", "--version"], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
	except FileNotFoundError:
		raise RuntimeError("sbatch command not found. Please make sure it is available in the PATH")

	return config_dict


def validate_pdb(filepath):
	'''
	Validate the PDB file by checking if it exists and contains ATOM records
	:param filepath: Path to the PDB file
	:return: True if the file is valid, otherwise raise an exception
	'''
	if not os.path.exists(filepath) or not os.path.isfile(filepath):
		raise RuntimeError(f"{filepath} is not a readable file")
	
	with open(filepath, 'r') as f:
		count_atom_lines = sum(1 for line in f if line.startswith("ATOM"))
	
	if count_atom_lines > 0:
		return True
	else:
		raise RuntimeError(f"{filepath} is not a valid PDB file")


def reset_jobid(job_id):
	'''
	Return the job ID if it is not None, otherwise return -1.
	This is for the case when only certain steps are run.
	By submitting -1 as the dependency, the job will not wait for any other job.
	:param job_id: Job ID
	:return: Job ID if not None, otherwise -1
	'''

	return job_id if job_id else -1


def submit_job(script, args=[], dependency=None, slurm_opts=None):
	'''
	Submit a job to Slurm with optional dependencies.
	:param script: Path to the script to run
	:param args: List of arguments to pass to the script
	:param dependency: Job ID to depend on
	:param slurm_opts: Dictionary of additional Slurm options
	'''

	cmd = ['sbatch']
	if dependency:
		cmd.append(f'--dependency={dependency}')
	if slurm_opts:
		cmd.extend(slurm_opts)
	cmd.append(script)
	cmd.extend(args) if len(args) > 0 else None
	result = subprocess.run(" ".join([str(item) for item in cmd]), capture_output=True, text=True, check=True, shell=True)

	return result.stdout.strip().split()[-1]


def create_motif_list(rec_name, pdb_list_file):
	motif_files = glob.glob(rf"???_{rec_name}.pdb")
	with open(pdb_list_file, "w") as f:
		f.writelines("\n".join(motif_files))
		# add a new line
		f.write("\n")


# Function to run createPDS directly
def run_createPDS(pdb_list_file, output=None):
	'''
	Run createPDS with the provided arguments
	:param master_dir: Path to the MASTER directory with the executables
	:param pdb_list_file: Path to the file containing the list of PDB files
	:param output: Path to the output file or None to print to stdout
	return: None
	'''
	cmd = [f"createPDS", "--type", "query", "--pdbList", pdb_list_file]
	if output:  # If you want to redirect output to a file or /dev/null
		with open(output, 'w') as f:
			subprocess.run(' '.join(cmd), stdout=f, check=True, stderr=f, shell=True)
	else:
		subprocess.run(' '.join(cmd), check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL, shell=True)


def parse_args():
	'''
	Parse command line arguments and return the parsed arguments
	return: Namespace object containing the parsed arguments
	'''
	parser = argparse.ArgumentParser(
		description="PatchMAN performs search on existing monomers and complexes with structural motifs extracted from the query receptor and extract complementary fragments to be used as templates for peptide-protein interactions.")
	parser.add_argument("receptor", type=str, help="PDB file with the receptor")
	parser.add_argument("peptide_sequence", type=str,
	                    help='Peptide sequence can include modified residues in "GFK[SER:phosphorylated]RAD" format.')
	parser.add_argument("-m", "--min_rec_bb", action="store_true",
	                    help="Minimize receptor backbone (default: false)")
	parser.add_argument("-g", "--log_file", type=str, default="stdout", help="Log file (Default is stdout)")
	parser.add_argument("-e", "--error_log_file", type=str, default="stderr", help="Error log file (Default is stderr)")
	parser.add_argument("-n", "--job_name", type=str, default="PatchMAN_JOB", help="Job name (Default: PatchMAN_JOB)")
	parser.add_argument("-v", "--verbose", action="store_true", help="Print information about the job")
	parser.add_argument("-w", "--work_dir", type=Path, default=Path.cwd(),
	                    help="Working directory (Default: current directory)")
	parser.add_argument("-a", "--native_pdb", type=Path,
	                    help="Native PDB if exists. Needs to be exactly the same chains and lengths for both receptor and peptide")
	parser.add_argument("-t", "--nstruct", type=int, default=1,
	                    help="Number of structures to generate (Default: 1)")
	parser.add_argument("-c", "--master_cutoff", type=float, default=1.5, help="Master cutoff (Default: 1.5)")
	parser.add_argument("-o", "--hotspot_mode", action="store_true",
	                    help="Hotspot mode, only the residues in focus will be used as patch centers (Default: false)")
	parser.add_argument("-s", "--mask_pdb", type=Path,
	                    help="PDB file with residues that should or should not be in the binding site (type: pdb file, Default: None)")
	parser.add_argument("-l", "--list_of_residues", type=str,
	                    help='List of residues that should or should not be in the binding site (type: string, Default: None, format: "A23,A12")')
	parser.add_argument("-f", "--focus_mode", action="store_true",
	                    help="Focus mode - the residues should from -l or -s should form the binding site (type: boolean, Default: False, masking mode)")
	parser.add_argument("-p", "--steps", type=str,
	                    help="Steps to run between (Default: 1-6, 1: split to motifs, 2: prepack receptor, 3: run MASTER, 4: extract templates,  5: FlexPepDock, 6: clustering and finalizing)")
	parser.add_argument("-b", "--benchmark", type=str, default=None,
	                    help="Benchmark mode, use the benchmark mode of FlexPepDock with the file provided (Default: None)")
	parser.add_argument("-u", "--cluster_radius", type=str, default="2",
	                    help="Cluster radius for clustering with Rosetta (Default: None)")
	
	args = parser.parse_args()
	
	# also validate working directory and create if it doesn't exist
	try:
		if not os.path.exists(args.work_dir):
			os.makedirs(args.work_dir)
	except Exception as e:
		raise ValueError(f"Working directory '{args.work_dir}' is not accessible and cannot be created: {e}")
	# add full path to the environment variable
	os.environ['work_dir'] = str(args.work_dir.resolve())
	
	# Additional processing for '-p' argument to split it into step_from and step_to
	step_from, step_to = 1, 6
	if args.steps:
		step_range = args.steps.split("-")
		step_from = int(step_range[0])
		step_to = int(step_range[1]) if step_range[1] != '' else 6
		
	args.steps_range = range(step_from, step_to + 1)
	
	# Validate inputs before proceeding with processing
	if not os.path.exists(args.receptor):
		raise ValueError(f"Receptor file '{args.receptor}' does not exist")
	if not os.path.isfile(args.receptor):
		raise ValueError(f"Receptor '{args.receptor}' is not a readable file")
	else:
		args.receptor = Path(args.receptor).resolve()
		# Validate the PDB file
		validate_pdb(args.receptor)
	
	# Ensure the receptor and native pdb are paths
	if args.native_pdb:
		args.native_pdb = args.native_pdb.resolve()
	if args.mask_pdb:
		args.mask_pdb = args.mask_pdb.resolve()

	# Additional arguments for FlexPepDock - store them in args for simplicity
	args.fpd_args = []
	if args.min_rec_bb:
		args.fpd_args.append(f"-m")
	if args.nstruct:
		args.fpd_args.append(f"-t {args.nstruct}")
	if args.native_pdb:
		args.fpd_args.append(f"-a {args.native_pdb}")
	return args


def prepare_pdb(filepath):
	'''
	Copy the PDB file to the current working directory and return the new path
	:param filepath: Path to the PDB file
	:return: New path to the copied PDB file
	'''
	validate_pdb(filepath)
	new_file = os.path.join(os.getcwd(), os.path.basename(filepath))
	shutil.copy(filepath, new_file)
	return new_file


def prepare_and_set_filenames(receptor_path):
	"""
	Prepares the receptor PDB file, extracts the base name, and sets up additional PDB filenames.

	Args:
	receptor_path: str or Path, path to the receptor PDB file.

	Returns:
	A tuple containing paths/names for the prepared receptor, clean receptor, and ppk receptor.
	"""
	
	# Prepare receptor PDB file
	receptor = prepare_pdb(receptor_path)
	receptor_base = os.path.basename(receptor)
	
	# Extract the base name without extension
	rec_name = os.path.splitext(receptor_base)[0]
	clean_rec = f"{rec_name}.clean.pdb"
	ppkrec = f"{rec_name}.clean.ppk.pdb"
	
	return receptor, receptor_base, clean_rec, ppkrec


def prepare_mask_and_focus(args):
	"""
	Prepares the mask and focus based on the provided inputs. Raises an error if both focus list and mask focus are provided.

	Args:
	focus_list: str or None, list of residues for focus.
	mask_focus: str or Path, path to the PDB file for mask focus.

	Returns:
	focus_mask_args: str, command-line arguments for focus and mask.

	Raises:
	SystemExit: If both focus list and mask focus are provided.
	"""
	focus_mask_args = ""
	
	focus_list = args.list_of_residues
	mask_focus = args.mask_pdb
	
	if focus_list and mask_focus:
		print("Provide either PDB file or list for mask/focus residues!", file=sys.stderr)
		sys.exit(1)
	elif mask_focus:
		mask_focus_prepared = prepare_pdb(mask_focus)  # Assumes prepare_pdb returns a str or Path
		focus_mask_args = f"-s {mask_focus_prepared}"
	
	return focus_mask_args


def validate_peptide_sequence(peptide_sequence):
    """
    Validates the given peptide sequence.

    Args:
    peptide_sequence: str, the peptide sequence to validate.

    Returns:
    The validated peptide sequence.

    Raises:
    ValueError: If the peptide sequence is invalid.
    """
    pep_sequence_to_validate = re.sub(r'\[[A-Z]{3,4}:[a-z]+\]', '', peptide_sequence)
    if not re.match(r"^[ARNDCEQGHILKMFPSTWYV]+$", pep_sequence_to_validate):
        raise ValueError(f"Not a peptide sequence: '{peptide_sequence}'")
    return peptide_sequence



def prepare_benchmark_mode(args, receptor_base):
	"""
	Prepare for benchmark mode by setting up the benchmark file and returning the path.

	Args:
	args: Namespace, parsed command-line arguments.

	Returns:
	benchmark_file: str, path to the benchmark file.
	"""
	if args.benchmark:
		# Read the benchmark file and extract relevant lines
		with open(args.benchmark, 'r') as file:
			for line in file:
				# Search for the first 4 letters of the receptor_base case-insensitively
				if re.search(receptor_base[:4], line, re.IGNORECASE):
					# Split the line by '|' and write to benchmark_file
					with open('benchmark_file', 'w') as outfile:
						outfile.write(line.replace('|', '\n'))
					break  # Assuming only the first match is needed, similar to grep -m 1
	
		# Filter out entries from the db_list_30nonred that are listed in benchmark_file
		with open(f'{os.environ["PROTOCOL_ROOT"]}/db_list_30nonred', 'r') as db_file, \
				open('benchmark_file', 'r') as benchmark_file, \
				open('custom_db_list_30nonred', 'w') as output_file:
			benchmark_entries = set(benchmark_file.read().splitlines())
			for db_line in db_file:
				# Check if any line from db_list_30nonred is not in the benchmark entries
				if db_line.strip() not in benchmark_entries:
					output_file.write(db_line)
		
		# Adjust master_args and fpd_args
		args.master_args = 'db_list_30nonred'
		args.fpd_args.append(f"{fpd_args} -b")
	else:
		args.master_args = f'{os.environ["PROTOCOL_ROOT"]}/db_list_30nonred'
	return args