#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=unmapped_%J_stdout.txt
#SBATCH --error=unmapped_%J_stderr.txt
#SBATCH --job-name=unmapped
# 

bash /home/biol726306/project/scripts/sort_bam/unmapped.sh