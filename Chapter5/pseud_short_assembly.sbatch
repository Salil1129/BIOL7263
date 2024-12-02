#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem 80G
#SBATCH --output=pseud_short_%J_stdout.txt
#SBATCH --error=pseud_short_%J_stderr.txt
#SBATCH --job-name=pseud_short
# 

bash /home/biol726306/chapter5/pseud_short_assembly.sh