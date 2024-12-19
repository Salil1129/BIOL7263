#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=tags_%J_stdout.txt
#SBATCH --error=tags_%J_stderr.txt
#SBATCH --job-name=add_tags
# 

bash /home/biol726306/project/scripts/sort_bam/add_tags.sh