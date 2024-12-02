# map
bwa mem -t 6 /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/hybrid_assembly \
/scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_1.fq.gz \
/scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/SRR491287_trimmed_reads_val_2.fq.gz \
-o /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_illumina.sam

# convert to bam
samtools view -bS /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseud_illumina.sam > /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina.bam

# sort
samtools sort -o /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina.bam

# index
samtools index /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam

# stats
samtools flagstat /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.bam > /scratch/biol726306/BIOL7263_Genomics/pseudomonas_gm41/mapping_to_assembly/pseudo_illumina_sorted.stats