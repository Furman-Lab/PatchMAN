#!/usr/bin/env python3

import os
import subprocess
import datetime

# Define paths and commands
python_cmd = os.environ.get("PYTHON", "python3")  # Default to "python3" if not found
protocol_root = os.environ.get("PROTOCOL_ROOT", "/path/to/protocol/root")  # Replace with default path if needed
rosetta_bin = os.environ.get("ROSETTA_BIN", "/path/to/rosetta/bin")  # Replace with default path if needed

# Sort by reweighted_sc
subprocess.run([python_cmd, os.path.join(protocol_root, "bin/printScoreFile_byHeader.py"),
				"score.sc", "reweighted_sc", "I_sc", "rmsBB", "rmsBB_if", "description"], stdout=subprocess.PIPE)
subprocess.run(["gawk", "{print $2, $3, $4, $5, $6}", "short.sc"], stdout=subprocess.PIPE)
subprocess.run(["sort", "-nk1", "short.sc"], stdout=subprocess.PIPE)

# Filter possible disulf bridges
tmp_file = "tmp"
with open("sorted.sc", "r") as sorted_file, open(tmp_file, "w") as tmp_output:
	for line in sorted_file:
		fields = line.split()
		score = fields[0].split('.')[0]
		if int(score) < -1000000:
			tmp_output.write(fields[-1] + "\n")

if os.path.exists(tmp_file):
	all_filter = "|".join(line.strip() for line in open(tmp_file))
	subprocess.run(["grep", "-v", "-E", all_filter, "sorted.sc"], stdout=open("tmp2", "w"))
	os.rename("tmp2", "sorted.sc")
	os.remove(tmp_file)

# Take top X percent structures
num_str = 0
with open("sorted.sc", "r") as sorted_file:
	nlines = sum(1 for _ in sorted_file)
	if nlines > 1000:
		num_str = nlines // 100
	elif nlines > 100:
		num_str = nlines // 10
	else:
		num_str = nlines // 2

subprocess.run(["head", f"-{num_str}", "sorted.sc"], stdout=open("top1percent", "w"))

# Run clustering
os.makedirs("clustering", exist_ok=True)
os.chdir("clustering")

pdb = next(iter(os.listdir("..")), None)
len_cmd = subprocess.run(["grep", "CA", pdb], stdout=subprocess.PIPE)
len_result = len_cmd.stdout.decode().strip()
len_val = int(len_result) if len_result.isdigit() else 0
peps = sys.argv[1]
plen = len(peps)
radius_arg = sys.argv[2]
actualR_cmd = subprocess.run(["date", "+%s"], stdout=subprocess.PIPE)
actualR_val = int(actualR_cmd.stdout.decode().strip()) * (plen / len_val) * float(radius_arg)
print(f"Actual radius is {actualR_val}")

# Calling mpiserialization to only have to compile Rosetta with these 2 extras
subprocess.run([os.path.join(rosetta_bin, "cluster.mpiserialization.linuxgccrelease"),
					"-in:file:silent", "../decoys.silent", "-in:file:silent_struct_type", "binary",
					"-cluster:radius", str(actualR_val), "-in:file:fullatom",
					"-tags", *open("../top1percent").read().split()], stdout=open("clog", "w"))

