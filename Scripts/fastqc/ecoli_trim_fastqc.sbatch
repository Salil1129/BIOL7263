#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=ecoli_fastqc_%J_stdout.txt
#SBATCH --error=ecoli_fastqc_%J_stderr.txt
#SBATCH --job-name=ecoli_fastqc
# 

bash /home/biol726306/BIOL7263_Genomics/scripts/fastqc/ecoli_trim_fastqc.sh