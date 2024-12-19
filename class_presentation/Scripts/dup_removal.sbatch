#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=dup_removal_%J_stdout.txt
#SBATCH --error=dup_removal_%J_stderr.txt
#SBATCH --job-name=dup_removal
# 

bash /home/biol726306/project/scripts/sort_bam/dup_removal.sh