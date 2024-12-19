#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 20G
#SBATCH --output=project_trim_%J_stdout.txt
#SBATCH --error=project_trim_%J_stderr.txt
#SBATCH --job-name=project_trim
# 

bash /home/biol726306/project/scripts/trim_galore/project_trim.sh