#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=sort_read_name_%J_stdout.txt
#SBATCH --error=sort_read_name_%J_stderr.txt
#SBATCH --job-name=sort_bam_read
# 

bash /home/biol726306/project/scripts/sort_bam/sort_read_name.sh