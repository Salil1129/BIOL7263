#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 100G
#SBATCH --output=pseud_long_%J_stdout.txt
#SBATCH --error=pseud_long_%J_stderr.txt
#SBATCH --job-name=pseud_long
# 

bash /home/biol726306/chapter5/pseudo_long_assembly.sh