#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 20G
#SBATCH --output=sort_chr_%J_stdout.txt
#SBATCH --error=sort_chr_%J_stderr.txt
#SBATCH --job-name=sort_chr_ag
# 

bash /home/biol726306/project/scripts/sort_bam/sort_chr_ag.sh