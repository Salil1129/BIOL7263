#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=project_fastqc_%J_stdout.txt
#SBATCH --error=project_fastqc_%J_stderr.txt
#SBATCH --job-name=project_fastqc
# 

bash /home/biol726306/project/scripts/fastqc/project_fastqc.sh