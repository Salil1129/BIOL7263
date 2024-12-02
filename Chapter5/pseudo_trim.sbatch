#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 80G
#SBATCH --output=pseudo_trim_%J_stdout.txt
#SBATCH --error=pseudo_trim_%J_stderr.txt
#SBATCH --job-name=pseudo_trim
# 

bash /home/biol726306/chapter5/pseudo_trim.sh