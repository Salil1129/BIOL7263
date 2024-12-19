#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 80G
#SBATCH --output=spades_rna_%J_stdout.txt
#SBATCH --error=spades_rna_%J_stderr.txt
#SBATCH --job-name=spades_rna
# 

bash /home/biol726306/project/scripts/sort_bam/spades_rna.sh