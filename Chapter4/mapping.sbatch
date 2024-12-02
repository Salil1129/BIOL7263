#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=reads_mapping_%J_stdout.txt
#SBATCH --error=reads_mapping_%J_stderr.txt
#SBATCH --job-name=reads_mapping
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/mapping.sh