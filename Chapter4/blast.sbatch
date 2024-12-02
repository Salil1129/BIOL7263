#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem 16G
#SBATCH --output=blast_%J_stdout.txt
#SBATCH --error=blast_%J_stderr.txt
#SBATCH --job-name=blast
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/blast.sh