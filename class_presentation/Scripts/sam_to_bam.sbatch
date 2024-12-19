#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 64G
#SBATCH --output=bam_to_sam_%j.out.txt
#SBATCH --error=bam_to_sam_%j.err.txt
#SBATCH --job-name=sam_to_bam
# 

bash /home/biol726306/project/scripts/sam_to_bam/sam_to_bam.sh