#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=bam_to_fasta_%J_stdout.txt
#SBATCH --error=bam_to_fasta_%J_stderr.txt
#SBATCH --job-name=bam_to_fasta
# 

bash /home/biol726306/project/scripts/sort_bam/bam_to_fasta.sh