#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=unmapped_fastqc_%J_stdout.txt
#SBATCH --error=unmapped_fastqc_%J_stderr.txt
#SBATCH --job-name=unmapped_fastqc
# 

bash /home/biol726306/project/scripts/sort_bam/unmapped_fastqc.sh