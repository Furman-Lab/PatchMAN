#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=13:00:00
#SBATCH --mem=1G
#SBATCH --module="singularity"

${PYTHON_SINGULARITY}/prepack_receptor.py $1
