#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=2
#SBATCH --mem 16G
#SBATCH --output=qualimap_%J_stdout.txt
#SBATCH --error=qualimap_%J_stderr.txt
#SBATCH --job-name=qualimap
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/qualimap.sh