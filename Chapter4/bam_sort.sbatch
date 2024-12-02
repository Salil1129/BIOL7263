#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=bam_sort_%J_stdout.txt
#SBATCH --error=bam_sort_%J_stderr.txt
#SBATCH --job-name=bam_sort
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/bam_sort.sh