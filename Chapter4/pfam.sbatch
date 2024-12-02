#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem 16G
#SBATCH --output=pfam_%J_stdout.txt
#SBATCH --error=pfam_%J_stderr.txt
#SBATCH --job-name=pfam
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/pfam.sh