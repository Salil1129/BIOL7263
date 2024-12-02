#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 80G
#SBATCH --output=illumina_map_%J_stdout.txt
#SBATCH --error=illumina_map_%J_stderr.txt
#SBATCH --job-name=illumina_map
# 

bash /home/biol726306/chapter5/illumina_map.sh