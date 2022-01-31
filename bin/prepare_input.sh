#!/bin/bash
#SBATCH --job-name=master
#SBATCH --time=13:00:00
#SBATCH --mem=1G

${PYTHON_SINGULARITY}/prepack_receptor.py $1
