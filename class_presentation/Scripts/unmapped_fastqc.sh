mkdir -p /scratch/biol726306/Project/fastqc_output/unmapped_fastqc

fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc
fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc