diamond blastx --threads 12 --outfmt 6 qseqid sseqid length pident evalue stitle -k 4 -d /home/biol726306/diamond\ database/viral_proteins.dmnd -q /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta -o /scratch/biol726306/Project/diamond_blastx/transcripts_blastx.tsv