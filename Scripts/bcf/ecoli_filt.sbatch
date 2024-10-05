#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_call_%J_stdout.txt
#SBATCH --error=ecoli_call_%J_stderr.txt
#SBATCH --job-name=ecoli_call
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/bcf/ecoli_filt.sh