#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=sort_bam_%J_stdout.txt
#SBATCH --error=sort_bam_%J_stderr.txt
#SBATCH --job-name=sort_bam
# 

bash /home/biol726306/project/scripts/sort_bam/sort_bam.sh