#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=sam_to_bam_%J_stdout.txt
#SBATCH --error=sam_to_bam_%J_stderr.txt
#SBATCH --job-name=sam_to_bam
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/assembly/sam_tobam.sh