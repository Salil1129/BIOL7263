#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=reads_index_%J_stdout.txt
#SBATCH --error=reads_index_%J_stderr.txt
#SBATCH --job-name=reads_index
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/reads_index.sh