#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 1024
#SBATCH --output=test_%J_stdout.txt
#SBATCH --error=test_%J_stdout.txt
#SBATCH --job-name=test
# 

bash test.sh