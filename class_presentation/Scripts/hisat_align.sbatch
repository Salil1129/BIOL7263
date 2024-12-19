#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 64G
#SBATCH --output=hisat2_%j.out.txt
#SBATCH --error=hisat2_%j.err.txt
#SBATCH --job-name=hisat2_align
# 

bash /home/biol726306/project/scripts/hisat2/hisat_align.sh