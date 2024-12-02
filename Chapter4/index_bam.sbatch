#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=index_bam_%J_stdout.txt
#SBATCH --error=index_bam_%J_stderr.txt
#SBATCH --job-name=index_bam
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/index_bam.sh