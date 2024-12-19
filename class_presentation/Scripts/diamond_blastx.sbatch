#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 20G
#SBATCH --output=blastp_diamond_%J_stdout.txt
#SBATCH --error=blastp_diamond_%J_stderr.txt
#SBATCH --job-name=blastp
# 

bash /home/biol726306/project/scripts/diamond_blastx/diamond_blastx.sh
