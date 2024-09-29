mkdir -p /scratch/biol726306/BIOL7263_Genomics/trimmedfastqc_output

fastqc /scratch/biol726306/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_1.fq.gz -o /scratch/biol726306/BIOL7263_Genomics/trimmedfastqc_output
fastqc /scratch/biol726306/BIOL7263_Genomics/sequencing_data/ecoli/trimmed_reads_val_2.fq.gz -o /scratch/biol726306/BIOL7263_Genomics/trimmedfastqc_output