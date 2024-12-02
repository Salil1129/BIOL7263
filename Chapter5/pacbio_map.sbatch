#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 80G
#SBATCH --output=pacbio_map_%J_stdout.txt
#SBATCH --error=pacbio_map_%J_stderr.txt
#SBATCH --job-name=pacbio_map
# 

bash /home/biol726306/chapter5/pacbio_map.sh