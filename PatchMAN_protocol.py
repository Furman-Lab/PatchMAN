#!/usr/bin/env python3

import os
import argparse
import subprocess
import shutil
import re
import datetime
from bin.split_to_motifs import main


def die(message):
	print(f"\nERROR: {message}\n", file=sys.stderr)
	sys.exit(1)


def validate_pdb(filepath):
	if not os.path.exists(filepath) or not os.path.isfile(filepath):
		die(f"{filepath} is not a readable file")
	
	with open(filepath, 'r') as f:
		count_atom_lines = sum(1 for line in f if line.startswith("ATOM"))
	
	if count_atom_lines > 0:
		return True
	else:
		die(f"{filepath} is not a valid PDB file")


def zero_jobid(job_id):
	return job_id if job_id else -1


def prepare_pdb(filepath):
	validate_pdb(filepath)
	new_file = os.path.join(os.getcwd(), os.path.basename(filepath))
	shutil.copy(filepath, new_file)
	return new_file


def print_jobid(job_type, job_id):
	if verbose:
		print(f"DEBUG| {job_type} JOBID: {job_id}")


# Parse command-line arguments
parser = argparse.ArgumentParser(description="The PatchMAN protocol")
parser.add_argument("receptor", help="PDB file with the receptor")
parser.add_argument("peptide_sequence", help="Peptide sequence")
parser.add_argument("-g", "--logs_dir", help="Log file directory")
parser.add_argument("-c", "--master_cutoff", type=float, default=1.5, help="Master cutoff")
parser.add_argument("-w", "--work_dir", default=os.getcwd(), help="Working directory")
parser.add_argument("-t", "--nstruct", type=int, help="Number of structures to generate")
parser.add_argument("-m", "--min_rec_bb", action="store_true", help="Minimize receptor backbone")
parser.add_argument("-p", "--step", type=int, default=1, help="Step to start from")
parser.add_argument("-v", "--verbose", action="store_true", help="Print information about the job")
parser.add_argument("-s", "--mask", help="Mask file with residues not in the binding site")
parser.add_argument("-f", "--focus", help="Focus mask with residues in the binding site")
# Add more argument definitions here...

args = parser.parse_args()

# Create output directory if it does not exist and change to it
os.makedirs(args.work_dir, exist_ok=True)
os.chdir(args.work_dir)

# VALIDATE INPUT
if not os.path.exists(args.receptor):
	die(f"Receptor file '{args.receptor}' does not exist")
if not os.path.isfile(args.receptor):
	die(f"Receptor '{args.receptor}' is not a readable file")

if not os.path.exists(args.work_dir):
	die(f"Working directory '{args.work_dir}' does not exist")
if not os.path.isdir(args.work_dir):
	die(f"'{args.work_dir}' is not a directory")

# Prepare receptor
receptor = prepare_pdb(args.receptor)
receptor_base = os.path.basename(receptor)

# Prepare input peptide
pep_sequence_to_validate = re.sub(r'\[[A-Z]{3,4}:[a-z]+\]', '', args.peptide_sequence)
if not re.match(r"^[ARNDCEQGHILKMFPSTWYV]+$", pep_sequence_to_validate):
	die(f"Not a peptide sequence: '{args.peptide_sequence}'")
pep_sequence = args.peptide_sequence

# Prepare mask and focus if provided
mask = args.mask if args.mask else ""
focus = args.focus if args.focus else ""

if mask and focus:
	die("Cannot provide both mask and focus")

if mask:
	mask = prepare_pdb(mask)
	mask_option = "-m " + mask
else:
	mask_option = ""

if focus:
	focus = prepare_pdb(focus)
	focus_option = "-f " + focus
else:
	focus_option = ""

# Prepare job
job_name = re.sub(r'[\t\n ]+', '_', args.job_name)
rec_name = os.path.splitext(os.path.basename(receptor))[0]
clean_rec = rec_name + ".clean.pdb"
ppkrec = rec_name + ".clean.ppk.pdb"

# Step 1: Split to motifs
if args.step <= 1:
	# Prepare PDB motifs and create MASTER pds files
	os.system(
		f"python {args.bin_dir}/split_to_motifs.py -i {receptor} {'-v' if args.verbose else ''} {'-m ' + mask if args.mask else ''} {'-f ' + focus if args.focus else ''}")
	
	motif_files = [f for f in os.listdir() if re.match(r'^[A-Z0-9_]{3}_' + receptor_base + r'$', f)]
	with open("motif_list", "w") as f:
		f.write('\n'.join(motif_files))
	
	os.system(f"{args.master_dir}/createPDS --type query --pdbList motif_list > /dev/null")

# Step 2: Prepack receptor
if args.step <= 2:
	prepack_receptor_jid = subprocess.run(
		["sbatch", "--job-name=prepack_receptor", "--get-user-env", "--time=90:00:00", "--mem=1G", "--module=singularity",
		 "--mem=1600m", f"--wrap='{args.bin_dir}/prepack_receptor.py {clean_rec}'"],
		stdout=subprocess.PIPE, text=True).stdout.split()[-1]

# Step 3: Run MASTER
if args.step <= 3:
	prepack_receptor_jid = zero_jobid(prepack_receptor_jid)
	print_jobid("PREPACK", prepack_receptor_jid)
	
	run_master_jid = subprocess.run(
		["sbatch", f"--dependency=afterok:{prepack_receptor_jid}",
		 f"--array=0-1000%50", f"{args.bin_dir}/run_master.sh {str(args.master_cutoff)}"],
		stdout=subprocess.PIPE, text=True).stdout.split()[-1]

# Step 4: Extract templates
if args.step <= 4:
	run_master_jid = zero_jobid(run_master_jid)
	print_jobid("MASTER", run_master_jid)
	
	extract_templates_jid = subprocess.run(
		["sbatch", "--array=0-1000%50", "--dependency=afterany:" + run_master_jid,
		 f"{args.bin_dir}/run_extract_templates.sh {pep_sequence}, {ppkrec}"],
		stdout=subprocess.PIPE, text=True).stdout.split()[-1]

# Step 5: FPD
if args.step <= 5:
	extract_templates_jid = zero_jobid(extract_templates_jid)
	print_jobid("EXTRACT TEMPLATES", extract_templates_jid)
	
	fpd_args = []
	if args.nstruct is not None:
		fpd_args.append("-t")
		fpd_args.append(str(args.nstruct))
	
	subprocess.run(
		["sbatch", f"--dependency=afterany:{str(extract_templates_jid)}",
		 f"--chdir={os.getcwd()}", f"--job-name=fpd {args.add_sbatch}",
		 f"{args.bin_dir}/fpd.sh -u {clean_rec} -m {str(args.min_rec_bb)} {fpd_args}"])

# Step 6: Clustering & Finalizing
if args.step <= 6:
	fpd_jid = zero_jobid(fpd_jid)
	print_jobid("FPD", fpd_jid)

	print("Running steps: clustering and finalizing")
	clustering_jid = subprocess.run(
		["sbatch", "--job-name=clustering", "--nice=8000", f"--chdir={os.getcwd()}",
		f"--dependency=aftercorr:{str(fpd_jid)} --kill-on-invalid-dep=yes",
		f"--get-user-env --wrap='{args.bin_dir}/cluster.py {pep_sequence} {str(args.cluster_radius)}'"],
		stdout=subprocess.PIPE, text=True).stdout.split()[-1]

	clustering_jid = zero_jobid(clustering_jid)
	subprocess.run(
		["sbatch", "--job-name=finalize", "--nice=8000", "--time=93:00:00",
		f"--mem=2000m --chdir={os.getcwd()}",
		f"--dependency=aftercorr:{str(clustering_jid)}", "--kill-on-invalid-dep=yes",
		f"--mem-per-cpu=2000", "--get-user-env", f"{args.bin_dir}/finalize.sh"])

# Print job information if verbose
if args.verbose:
	print("------------------------------------------------")
	print(datetime.datetime.now())
	print(f"Verbosity: {args.verbose}")
	print(f"Receptor is: {receptor}")
	print(f"Peptide sequence is: {pep_sequence}")
	print(f"Cluster radius: {args.cluster_radius}")
	print(f"Slurm job IDs: {run_master_jid} {extract_templates_jid} {fpd_jid} {clustering_jid} {finalize_jid}")
	print(f"Current working directory: {os.getcwd()}")
	print("------------------------------------------------")

print(f"SLURM_JID {finalize_jid}")