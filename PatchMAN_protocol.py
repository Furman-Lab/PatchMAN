#!/usr/bin/env python3
from bin.protocol_utils import *

# Gather all the arguments needed for the protocol
args = parse_args()
config = load_config('config.ini')
os.chdir(args.work_dir)

# Prepare inputs
pep_sequence = validate_peptide_sequence(args.peptide_sequence)
args.focus_mask_args = prepare_mask_and_focus(args)

# Prepare receptor and derived filenames
receptor, receptor_base, ppkrec = prepare_and_set_filenames(args.receptor)
args.fpd_args.append(f"-u {receptor}")

# Copy also native file into directory
if args.native_pdb:
    args.native_pdb = prepare_pdb(args.native_pdb)
    args.fpd_args.append(f"-a {args.native_pdb}")

# Prepare for benchmark mode if provided
args = prepare_benchmark_mode(args, receptor_base)

# Step 1: Split to motifs
pdb_list_file = "motif_list"
if 1 in args.steps_range:
    # Call split_to_motifs, with provided arguments
    # Execute Python script directly, not a Slurm job
    subprocess.run(['python3', f"{os.environ['BIN_DIR']}/split_to_motifs.py", "-i", receptor,
                    *args.focus_mask_args.split()], check=True)

# Calculate n_searches using motif_list
n_searches = sum(1 for line in open(pdb_list_file))
print(f"Number of searches: {n_searches}")

# Step 2 to Step 6: Rewriting using submit_job
# Example for Step 2 - adjust according to actual scripts and parameters needed
if 2 in args.steps_range:
    prepack_receptor_jid = submit_job(
        f"{os.environ['BIN_DIR']}/prepare_input.sh"
    )
else:
    prepack_receptor_jid = -1

# Step 3: Run MASTER
if 3 in args.steps_range:
    run_master_jid = submit_job(
        script=f"{os.environ['BIN_DIR']}/run_master.sh",
        args=[str(args.master_cutoff), args.master_args],
        dependency=f'afterok:{prepack_receptor_jid}',
        slurm_opts={f"--array=0-{n_searches}%50"}
    )
else:
    run_master_jid = -1

# Step 4: Extract templates
if 4 in args.steps_range:
    extract_templates_jid = submit_job(
        script=f"{os.environ['BIN_DIR']}/run_extract_templates.sh",
        args=[args.peptide_sequence, ppkrec, *args.focus_mask_args.split()],
        dependency=f'afterany:{run_master_jid}',
        slurm_opts={f"--array=0-{n_searches}%50"}
    )
else:
    extract_templates_jid = -1

# Step 5: FlexPepDock (FPD)
if 5 in args.steps_range:
    fpd_jid = submit_job(
        script=f"{os.environ['BIN_DIR']}/fpd.sh",
        args=args.fpd_args,
        dependency=f'afterany:{extract_templates_jid}',
        slurm_opts=[f"--chdir={os.getcwd()}", config['ADD_SBATCH']]
    )
else:
    fpd_jid = -1

# Step 6: Clustering & Finalizing
if 6 in args.steps_range:
    clustering_jid = submit_job(
        script=f"{os.environ['BIN_DIR']}/cluster.sh",
        args=[args.peptide_sequence, str(args.cluster_radius)],
        dependency=f'aftercorr:{fpd_jid}',
        slurm_opts={"--job-name=clustering", "--nice=8000", f"--chdir={os.getcwd()}"}
    )
else:
    clustering_jid = -1


# Print job IDs for debugging or tracking
print(f'Submitted receptor {receptor_base} with peptide: {args.peptide_sequence}')
print("Job IDs:", {"Prepack": prepack_receptor_jid,
                   "Run MASTER": run_master_jid,
                   "Extract Templates": extract_templates_jid,
                   "FPD": fpd_jid,
                   "Clustering & Finalizing": clustering_jid})
