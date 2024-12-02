#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem 16G
#SBATCH --output=blast_p_%J_stdout.txt
#SBATCH --error=blast_p_%J_stderr.txt
#SBATCH --job-name=blast_p
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/blast_p.sh