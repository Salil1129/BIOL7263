blastn -query "/scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta" \
       -db /scratch/biol726306/blast_database/virus_db \
       -out /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt \
       -outfmt "6 qseqid sseqid pident stitle" \
       -max_target_seqs 1 | sort -u > /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt