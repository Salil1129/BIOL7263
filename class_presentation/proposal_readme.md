**Project background**
Tomato (Solanum lycopersicum L.) is a night shade family crop with high economic value for its produce.
However, the tomato production is restricted by numerous diseases and, of all pathogens, viruses are known to significantly affect its production. 
Earlier virus detection techniques were traditional targeted detection techniques including ELISA, PCR, and Sanger sequencing which require prior information on viral genomes or information for serological properties of viral species. 
Presently, the lower cost of high-throughput sequencing (HTS) and the availability of the improved data analysis tools has made possible the study of a whole community of viruses (i.e., viromes), including unknown ones.
This advancement has escaped challenges of targeted detection of plant viruses and contributed useful ecological and epidemiological insights.

**Research question**
How to identify viral contigs using RNA-seq data derived from infected tomato samples?

**Objectives**
1.	Identification of viral contigs present in the RNA-seq data of infected tomato samples
2.	Detection of any novel or uncharacterized viruses if present in the RNA-seq data

**Data source**
The total RNA extracted from the tomato samples collected during the field survey was subjected to high-throughput sequencing (HTS) using the NextSeq 500/550 high-output kit v2.5 (Illumina, U.S.A.) at the Genomic Facility, Oklahoma State University (Stillwater, OK), as part of our lab's research

**Analysis**

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

6. To check the quality of reads, fastqc command was ran with:

- sbatch file:
[sbatch_file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/project_fastqc.sbatch)
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

[fastqc.sh](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/project_fastqc.sh)


```bach
mkdir -p /scratch/biol726306/Project/fastqc_output

fastqc /scratch/biol726306/Project/read_1.fastq.gz -o /scratch/biol726306/Project/fastqc_output
fastqc /scratch/biol726306/Project/read_2.fastq.gz -o /scratch/biol726306/Project/fastqc_output
```
- Per base sequence quality showed a drop in quality towards the 3’ end of reads.
- Per Sequence Quality Scores showed that the most of the reads are acceptable,
there may be a significant proportion of lower-quality reads in the dataset.
- Adapter contamination was also present in the data.

7. To trim the adapter sequences, trim galore was used:

- .sbatch file:
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/project_trim.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/project_trim.sh)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/hisat_align.sh)
```bach
hisat2 -x /home/biol726306/project/reference_genome/index_basename -1 /scratch/biol726306/Project/read_1.fastq.gz -2 /scratch/biol726306/Project/read_2.fastq.gz -S /scratch/biol726306/Project/project_output.sam
```

- sbatch file:
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/hisat_align.sbatch)
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=4
#SBATCH --mem 64G
#SBATCH --output=hisat2_%j.out.txt
#SBATCH --error=hisat2_%j.err.txt
#SBATCH --job-name=hisat2_align
# 

bash /home/biol726306/project/scripts/hisat2/hisat_align.sh
```

14. Conversion of SAM to BAM file:

- sh file
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sam_to_bam.sh)
```bach
samtools view -bS /scratch/biol726306/Project/hisat2_reads_to_reference/project_output.sam -o /scratch/biol726306/Project/sam_to_bam/project_output.bam
```
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sam_to_bam.sbatch)
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
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_bam.sbatch)
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=sort_bam_%J_stdout.txt
#SBATCH --error=sort_bam_%J_stderr.txt
#SBATCH --job-name=sort_bam
# 

bash /home/biol726306/project/scripts/sort_bam/sort_bam.sh
```

- sh file
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_bam.sh)
```bach
samtools sort /scratch/biol726306/Project/sam_to_bam/project_output.bam \
-o /scratch/biol726306/Project/sort_bam/project_output_sorted.bam
```
16. Sorting the BAM file by 'read name'

- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_read_name.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_read_name.sh)
```bach
samtools sort -n \
/scratch/biol726306/Project/sort_bam/project_output_sorted.bam \
-o /scratch/biol726306/Project/sort_bam/project_output_namesort.bam
```

17. Adding tags to file
- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/add_tags.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/add_tags.sh)
```bach
samtools fixmate -m /scratch/biol726306/Project/sort_bam/project_output_namesort.bam \
/scratch/biol726306/Project/sort_bam/project_output_namesort_fixmate.bam
```
18. Resorting by chromosomal co-ordinates
- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_chr_ag.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/sort_chr_ag.sh)
```bach
samtools sort -o /scratch/biol726306/Project/sort_bam/namesort_fixmate_sort_chr.bam \
/scratch/biol726306/Project/sort_bam/project_output_namesort_fixmate.bam
```

19. removal of duplicates
- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/dup_removal.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/dup_removal.sh)
```bach
samtools markdup -r /scratch/biol726306/Project/sort_bam/namesort_fixmate_sort_chr.bam \
/scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam 
```
20. Generate index file
- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/index_bam.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/index_bam.sh)
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
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/qualimap.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/qualimap.sh)
```bach
qualimap bamqc -outdir /scratch/biol726306/Project/sort_bam/bamqc \
-bam /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam \
-gff /home/biol726306/project/reference_genome/GCF_000188115.5_SL3.1_genomic.gff
```
23. Collection of unmapped reads

To look at the file
```bach
samtools view /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam | head -n 5
```

**To filter out unmapped reads**
- sh file
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/unmapped.sh)
```bach
samtools view -b -f 12 /scratch/biol726306/Project/sort_bam/dup_removal_project_markdup.bam -o /scratch/biol726306/Project/unmapped_assembly/unmapped.bam
```
-sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/unmapped.sbatch)
```bach
#!/bin/bash
#
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --mem 16G
#SBATCH --output=unmapped_%J_stdout.txt
#SBATCH --error=unmapped_%J_stderr.txt
#SBATCH --job-name=unmapped
# 

bash /home/biol726306/project/scripts/sort_bam/unmapped.sh
```
24. To convert .bam to .fasta files

- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/bam_to_fasta.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/bam_to_fasta.sh)
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
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/unmapped_fastqc.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/unmapped_fastqc.sh)
```bach
mkdir -p /scratch/biol726306/Project/fastqc_output/unmapped_fastqc

fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc
fastqc /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq -o /scratch/biol726306/Project/fastqc_output/unmapped_fastqc
```

26. Assembly using SPAdes RNA
- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/spades_rna.sbatch)

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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/spades_rna.sh)
```bach
spades.py --rna -t 20 -m 80 -o /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly -1 /scratch/biol726306/Project/unmapped_assembly/unmapped_r1.fastq -2 /scratch/biol726306/Project/unmapped_assembly/unmapped_r2.fastq
```
- This took ~45 minutes to run

27. Quality assessment of the assembly using quast tool
go to the directory containig "transcript.fasta" file and run this command:

```bach
quast.py --output-dir quast contigs.fasta
```
![quast report](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/quast_report.png)
28 BLAST search of the unmapped contigs using diamond database

- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/diamond_blastx.sbatch)

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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/diamond_blastx.sh)

```bach
diamond blastx --threads 12 --outfmt 6 qseqid sseqid length pident evalue stitle -k 4 -d /home/biol726306/diamond\ database/viral_proteins.dmnd -q /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta -o /scratch/biol726306/Project/diamond_blastx/transcripts_blastx.tsv
```

29. BLAST search of the unmapped contigs against nucleotide database

- sbatch file
[sbatch file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/blast_search.sbatch)
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
[sh file](https://github.com/Salil1129/BIOL7263/blob/main/class_presentation/Scripts/blast_search.sh)

```bach
blastn -query "/scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/transcripts.fasta" \
       -db /scratch/biol726306/blast_database/virus_db \
       -out /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt \
       -outfmt "6 qseqid sseqid pident stitle" \
       -max_target_seqs 1 | sort -u > /scratch/biol726306/Project/unmapped_assembly/spades_rna_assembly/blast_results.txt
```
***Conclusion***
This project has been both challenging and rewarding. 
I began by processing RNA-seq data from infected tomato sample, focusing on quality control and viral contig identification using different tools. 
The task of distinguishing viral sequences from the host genome was complicated, particularly in optimizing the scripts, but the discovery of Potato leafroll virus (PLRV), Tomato chlorosis virus (ToCV), and a potential novel virus made the effort worthwhile. 
This experience sharpened my RNA-seq analysis skills and deepened my understanding of virology, particularly in detecting and characterizing plant viruses. 
I successfully identified two known viruses, PLRV and ToCV, both significant pathogens affecting tomato crops, confirming that the infected samples contained known viral threats. 
More interestingly, I discovered a viral contig that could represent a novel plant virus, as it did not match any existing viral sequences in current databases, suggesting the possibility of a new virus potentially affecting tomato plants.
Further validation and research are needed to confirm whether this contig corresponds to an entirely new virus. 
Moving forward, future considerations include the validation of this novel virus, the characterization of the identified viruses, and reporting the findings.