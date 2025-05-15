# RNA-seq Tutorial 
This repository contains scripts for a full RNA-Seq analysis pipeline. Starting with raw sequencing reads, it undergoes quality control, read trimming, alignment, quantification, and statistical analysis to obtain differential expression results. 

## Pipeline Overview 

### 1. Quality Control on Raw Reads 
**FastQC** on raw FASTQ files  
**MultiQC** summary report 

### 2. Adapater Trimming 
**Trimmomatic** to trim adapter sequences 

### 3. Quality Control on Trimmed Reads 
**FastQC** on trimmed FASTQ files  
**MultiQC** summary report 

### 4. Generate STAR Genome Index 
**STAR** on downloaded reference genome 

### 5. Read Alignment 
**STAR** to align reads to reference genome 

### 6. Quality Control of Aligned Files 
**FastQC** on *BAM* files  
**MultiQC** summary report 

### 7. Quantification 
**FeatureCounts** to quantify gene expression 

### 8. Differential Expression Analysis 
**EdgeR** for statistical analysis  
Includes data normalization, exploratory analyses, and visualization techniques 

## 1. Quality Control on Raw Reads
**Inputs:** .fastq.gz files  
**Outputs:** .html & .zip QC reports; MultiQC summary report  
**Bash:**  
```
fastqc 0_raw_data/*.fastq.gz -o 1_qc_raw/  
multiqc 1_qc_raw/ -o 1_qc_raw/  
```

## 2. Adapter Trimming 
**Inputs:** _R1.fastq.gz & _R2.fastq.gz  
**Outputs:** Paired and unpaired trimmed *FASTQ* files  
**Bash:**  
trimmomatic PE -threads 8 -phred33 \  
  0_raw_data/sample_R1.fastq.gz 0_raw_data/sample_R2.fastq.gz \  
  2_trimmed_reads/sample_R1_paired.fastq.gz 2_trimmed_reads/sample_R1_unpaired.fastq.gz \  
  2_trimmed_reads/sample_R2_paired.fastq.gz 2_trimmed_reads/sample_R2_unpaired.fastq.gz \  
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \  
  SLIDINGWINDOW:4:20 \  
  MINLEN:36  

## 3. Quality Control on Trimmed Reads 
**Inputs:** _paired.fastq.gz  
**Outputs:** post-trimming QC reports; MultiQC summary  
**Bash:**  
fastqc 2_trimmed_reads/*_paired.fastq.gz -o 3_qc_trimmed/  
multiqc 3_qc_trimmed/ -o 3_qc_trimmed/  

## 4. Generate STAR Genome Index 
**Inputs:** Genome *FASTA* file; *GTF* annotation file  
**Outputs:** *STAR* genome index files  
**Bash:**  
STAR --runThreadN 8 \  
  --runMode genomeGenerate \  
  --genomeDir 4_star_index/ \  
  --genomeFastaFiles genome.fa \  
  --sjdbGTFfile annotation.gtf \  
  --sjdbOverhang 99  

## 5. Read Alignment 
**Inputs:** Trimmed paired *FASTQ* files, *STAR* genome index  
**Outputs:** *BAM* files; *STAR* log files  
**Bash:**  
STAR --runThreadN 8 \  
  --genomeDir 4_star_index/ \  
  --readFilesIn 2_trimmed_reads/sample_R1_paired.fastq.gz 2_trimmed_reads/sample_R2_paired.fastq.gz \  
  --readFilesCommand zcat \  
  --outFileNamePrefix 5_alignment/sample_ \  
  --outSAMtype BAM SortedByCoordinate  

## 6. Quality Control of Aligned Files  
**Inputs:** .bam files from STAR alignment  
**Outputs:** QC summary file  
**Bash:**  
fastqc 5_alignment/*.bam -o 6_qc_alignment/  
multiqc 6_qc_alignment/ -o 6_qc_alignment/  

## 7. Quantification  
**Inputs:** *BAM* files, *GTF* annotation files  
**Outputs:** Count matrix (.txt)  
**Bash:**  
featureCounts -T 8 -p \  
  -a annotation.gtf \  
  -o 7_featureCounts/featureCounts_matrix.txt \  
  5_alignment/*.bam  

## 8. Differential Expression Analysis  
**Inputs:** *FeatureCounts* count matrix  
**Outputs:**
**Bash:**
