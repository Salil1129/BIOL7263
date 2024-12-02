1. HTS data of a pooled tomato sample was downloaded into the scratch/biol726306/project directory using winscp (AR008_B)

2. The environment was activated using:
```bach
mamba activate /home/mbtoomey/.conda/envs/BIOL7263_Genomics
```

3. For the sake of simplicity, symbolic links were created using:
```bach
ln -s AR008-B_S186_R1_001.fastq.gz read_1.fastq.gz
ln -s AR008-B_S186_R2_001.fastq.gz read_2.fastq.gz
```
3.Reads were identified from the files as to work with sequence data require that the 'read 1' and 'read 2' files have the reads in the same order: However, this command was used first to identify the identifier: 
```bach
zcat read_1.fastq.gz | head -n 20
```
- Based on this, it was found out that @LH is the identifier in this data. Thus this command was used to identify order of the reads:
```bach
zcat read_1.fastq.gz | head | grep @LH
zcat read_2.fastq.gz | head | grep @LH
```
4. To get some more confidence, the above commands were repeated  using ‘tail’ instead of ‘head’ to compare reads at the end of the files.

```bach
zcat read_1.fastq.gz | tail | grep @LH
zcat read_2.fastq.gz | tail | grep @LH
```
5. To check if there is identical number of reads in each file, following command was used:
```bach
zcat read_1.fastq.gz | grep @LH | wc -l
zcat read_2.fastq.gz | grep @LH | wc -l
```
- It was found that both the files had identical number of reads viz. 25817149
6. fastqc command was ran with:

- sbatch file:
[sbach_file]("C:\Users\Dr Salil Jindal\OneDrive - University of Tulsa\Courses\Fall 2024\BIOL 7263\Proposal\project_fastqc.sbatch")
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=project_fastqc_%J_stdout.txt
#SBATCH --error=project_fastqc_%J_stderr.txt
#SBATCH --job-name=project_fastqc
# 

bash /home/biol726306/project/scripts/fastqc/project_fastqc.sh
```

- .sh file:

[fastqc.sh file]("C:\Users\Dr Salil Jindal\OneDrive - University of Tulsa\Courses\Fall 2024\BIOL 7263\Proposal\read_1_fastqc.html")


```bach
mkdir -p /scratch/biol726306/Project/fastqc_output

fastqc /scratch/biol726306/Project/read_1.fastq.gz -o /scratch/biol726306/Project/fastqc_output
fastqc /scratch/biol726306/Project/read_2.fastq.gz -o /scratch/biol726306/Project/fastqc_output
```

7. To trim the adapter sequences, trim galore was used:

- .sbatch file:
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 8G
#SBATCH --output=project_trim_%J_stdout.txt
#SBATCH --error=project_trim_%J_stderr.txt
#SBATCH --job-name=project_trim
# 

bash /home/biol726306/project/scripts/trim_galore/project_trim.sh
```
- sh file:
``` bach
trim_galore --paired --fastqc --gzip --cores 4 --length 100 /scratch/biol726306/Project/read_1.fastq.gz /scratch/biol726306/Project/read_2.fastq.gz --basename trimmed_reads -o /scratch/biol726306/Project/trim_galore/
```
8. After trimming, both the files were checked again for the identical number of reads using:
```bach
zcat trimmed_reads_val_1.fq.gz | wc -l
zcat trimmed_reads_val_2.fq.gz | wc -l
```
Both the files had identical number of reads viz. 87192560.
9. Downloading the reference genome using wget command:

- Open NCBI database homepage and look for the genome of interest which was tomato in this case.
- click on FTP and an index page will open.

[index_page](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/)

- Look for .fna.gz and .gff.gz files.

- For .fna.gz:
```bach
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.fna.gz
```

- For .gff.gz:
```bach
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/188/115/GCF_000188115.5_SL3.1/GCF_000188115.5_SL3.1_genomic.gff.gz

10. hisat2 was installed using 

```bach
conda create -n hisat2_env
```
and actiavted using:
```bach
conda activate hisat2_env
```
11. conversion from .fna.gz to .fna using
```bach
gunzip command
```

12. created index files using:
```bach
hisat2-build /home/biol726306/project/reference_genome/GCF_000188115.5_SL3.1_genomic.fna index_basename
```

13. Aligned index files of the reference genome to the reads:

- sh file:
```bach
hisat2 -x /home/biol726306/project/reference_genome/index_basename -1 /scratch/biol726306/Project/read_1.fastq.gz -2 /scratch/biol726306/Project/read_2.fastq.gz -S /scratch/biol726306/Project/project_output.sam
```

- sbatch file:
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 64G
#SBATCH --output=hisat2_%j.out.txt
#SBATCH --error=hisat2_%j.err.txt
#SBATCH --job-name=ecoli_vcf
# 

bash /home/biol726306/project/scripts/hisat2/hisat_align.sh
```

14. Conversion of SAM to BAM file:

- sh file
```bach
samtools view -bS /scratch/biol726306/Project/hisat2_reads_to_reference/project_output.sam -o /scratch/biol726306/Project/sam_to_bam/project_output.bam
```

- sbatch file
```bach
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
```
15. Soring the BAM file (by chromosomal co-ordinates)
- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=ecoli_view_%J_stdout.txt
#SBATCH --error=ecoli_view_%J_stderr.txt
#SBATCH --job-name=ecoli_view
# 

bash /home/biol726306/project/scripts/sort_bam/sort_bam.sh
```

- sh file
```bach
samtools sort /scratch/biol726306/Project/sam_to_bam/project_output.bam \
-o /scratch/biol726306/Project/sort_bam/project_output_sorted.bam
```
16. Sorting the BAM file by 'read name'

- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=sort_read_name_%J_stdout.txt
#SBATCH --error=sort_read_name_%J_stderr.txt
#SBATCH --job-name=sort_bam_read
# 

bash /home/biol726306/project/scripts/sort_bam/sort_read_name.sh
```
- sh file
```bach
samtools sort -n \
/scratch/biol726306/Project/sort_bam/project_output_sorted.bam \
-o /scratch/biol726306/Project/sort_bam/project_output_namesort.bam
```

17. Adding tags to file
- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=tags_%J_stdout.txt
#SBATCH --error=tags_%J_stderr.txt
#SBATCH --job-name=add_tags
# 

bash /home/biol726306/project/scripts/sort_bam/add_tags.sh
```
- sh file
```bach
samtools fixmate -m /scratch/biol726306/Project/sort_bam/project_output_namesort.bam \
/scratch/biol726306/Project/sort_bam/project_output_namesort_fixmate.bam
```
18. Resorting by chromosomal co-ordinates
- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 20G
#SBATCH --output=sort_chr_%J_stdout.txt
#SBATCH --error=sort_chr_%J_stderr.txt
#SBATCH --job-name=sort_chr_ag
# 

bash /home/biol726306/project/scripts/sort_bam/sort_chr_ag.sh
```

- sh file
```bach
samtools sort -o /scratch/biol726306/Project/sort_bam/namesort_fixmate_sort_chr.bam \
/scratch/biol726306/Project/sort_bam/project_output_namesort_fixmate.bam
```

19. removal of duplicates
- shbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=dup_removal_%J_stdout.txt
#SBATCH --error=dup_removal_%J_stderr.txt
#SBATCH --job-name=dup_removal
# 

bash /home/biol726306/project/scripts/sort_bam/dup_removal.sh
```

- sh file
```bach
samtools markdup -r /scratch/biol726306/Project/sort_bam/namesort_fixmate_sort_chr.bam \
/scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam 
```
20. Generate index file
- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=index_bam_%J_stdout.txt
#SBATCH --error=index_bam_%J_stderr.txt
#SBATCH --job-name=index_bam
# 

bash /home/biol726306/project/scripts/sort_bam/index_bam.sh
```

- sh file
```bach
samtools index /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam
```
- It stores the ouput file in the same folder as the input file but in a different format viz. ".bai"

21. Mapping statistics
```bach
samtools flagstat /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam
```
22. QualiMap

- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=qualimap_%J_stdout.txt
#SBATCH --error=qualimap_%J_stderr.txt
#SBATCH --job-name=qualimap
# 

bash /home/biol726306/project/scripts/sort_bam/qualimap.sh
```

- sh file
```bach
qualimap bamqc -outdir /scratch/biol726306/Project/sort_bam/bamqc \
-bam /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam \
-gff /home/biol726306/project/reference_genome/GCF_000188115.5_SL3.1_genomic.gff
```
23. Collection of unmapped reads
First look at the file:
```bach
samtools view /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam | head -n 5
```

**To filter out unmapped reads**
```bach
samtools view -b -f 12 /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam -o /scratch/biol726306/Project/unmapped_assembly/unmapped.bam
```
24. To convert .bam to .fasta files

- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=bam_to_fasta_%J_stdout.txt
#SBATCH --error=bam_to_fasta_%J_stderr.txt
#SBATCH --job-name=bam_to_fasta
# 

bash /home/biol726306/project/scripts/sort_bam/bam_to_fasta.sh
```

- sh file
```bach
bedtools bamtofastq -i /scratch/biol726306/Project/unmapped_assembly/unmapped.bam \
-fq /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq \
-fq2 /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq
```
- To check the number of reads
```bach
grep -c "^@LH" unmapped_r1.fastq unmapped_r2.fastq 
```
25. To check the quality
- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 8G
#SBATCH --output=unmapped_fastqc_%J_stdout.txt
#SBATCH --error=unmapped_fastqc_%J_stderr.txt
#SBATCH --job-name=unmapped_fastqc
# 

bash /home/biol726306/project/scripts/sort_bam/unmapped_fastqc.sh
```

- sh file
```bach
mkdir -p /scratch/biol726306/Project/fastqc_output/unmapped_fastqc

fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc
fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc
```

26. Assembly using SPAdes RNA
- batch file

```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 80G
#SBATCH --output=spades_rna_%J_stdout.txt
#SBATCH --error=spades_rna_%J_stderr.txt
#SBATCH --job-name=spades_rna
# 

bash /home/biol726306/project/scripts/sort_bam/spades_rna.sh
```

- sh file
```bach
spades.py --rna -t 20 -m 80 -o /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly -1 /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq -2 /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq
```
- This took ~45 minutes to run

27. Quality assessment of the assembly using quast tool
go to the directory containig "transcript.fasta" file and run this command:

```bach
quast.py --output-dir quast contigs.fasta
```
28 BLAST search of the unmapped contigs using diamond database

- sbatch file

```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem 20G
#SBATCH --output=blastp_diamond_%J_stdout.txt
#SBATCH --error=blastp_diamond_%J_stderr.txt
#SBATCH --job-name=blastp
# 

bash /home/biol726306/project/scripts/diamond_blastx/diamond_blastx.sh
```

- sh file

```bach
diamond blastx --threads 12 --outfmt 6 qseqid sseqid length pident evalue stitle -k 4 -d /home/biol726306/diamond\ database/viral_proteins.dmnd -q /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta -o /scratch/biol726306/Project/diamond_blastx/transcripts_blastx.tsv
```

29. BLAST search of the unmapped contigs against nucleotide database

- sbatch file
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 80G
#SBATCH --output=blast_search_%J_stdout.txt
#SBATCH --error=blast_search_%J_stderr.txt
#SBATCH --job-name=blast_search
# 

bash /home/biol726306/project/scripts/diamond_blastx/blast_search.sh
```
- sh file

```bach
blastn -query "/scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta" \
       -db /scratch/biol726306/blast_database/virus_db \
       -out /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt \
       -outfmt "6 qseqid sseqid pident stitle" \
       -max_target_seqs 1 | sort -u > /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt
```
