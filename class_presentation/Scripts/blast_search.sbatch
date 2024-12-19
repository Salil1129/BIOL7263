#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 80G
#SBATCH --output=blast_search_%J_stdout.txt
#SBATCH --error=blast_search_%J_stderr.txt
#SBATCH --job-name=blast_search
# 

bash /home/biol726306/project/scripts/diamond_blastx/blast_search.sh