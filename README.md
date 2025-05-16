# RNA-seq Tutorial 
This repository contains scripts for a full RNA-Seq analysis pipeline. Starting with raw sequencing reads, it undergoes quality control, read trimming, alignment, quantification, and statistical analysis to obtain differential expression results. 

## Pipeline Overview 

### 1. Quality Control on Raw Reads 
This steps checks the quality of the raw sequencing data to identify any issues such as low-quality scores or overrepresented sequences before processing. *MultiQC* generates a report of the aggregated QC results to help identify any issues in samples.  
**FastQC** on raw FASTQ files  
**MultiQC** summary report 

### 2. Adapater Trimming 
Here, any unwanted sequences that are added when sequencing such as adapter sequences or low-quality ends of reads are removed to improve alignment accuracy.  
**Trimmomatic** to trim adapter sequences 

### 3. Quality Control on Trimmed Reads 
Quality of reads are re-checked here to ensure that trimming actually improved read quality and did not introduce any new issues.  
**FastQC** on trimmed FASTQ files  
**MultiQC** summary report 

### 4. Generate STAR Genome Index 
In order to align reads, *STAR* needs a reference genome index so that it can locate what part of the genome the RNA reads came from.  
**STAR** on downloaded reference genome 

### 5. Read Alignment 
In this step, the reads are mapped to the reference genome, which provides the genes/regions that the RNA came from.  
**STAR** to align reads to reference genome 

### 6. Quantification 
Here, *FeatureCounts* counts how many reads mapped to each gene. The result is a matrix where each row is a gene and each column is a sample, with the data being the number of reads for each gene in each sample.  
**FeatureCounts** to quantify gene expression 

### 7. Quality Control of Aligned Files 
QC once again to see how well the reads mapped to genes and how many were counted. This further helps identify any disparities between or identify low-quality samples.  
**FastQC** on *BAM* files  
**MultiQC** summary report 

### 8. Differential Expression Analysis 
Here, statistical modeling is used to identify differentially expressed genes between experimental groups. Further, *EdgeR* can help visualize any trends between samples and genes.    
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
```
trimmomatic PE -threads 8 -phred33 \  
  0_raw_data/sample_R1.fastq.gz 0_raw_data/sample_R2.fastq.gz \  
  2_trimmed_reads/sample_R1_paired.fastq.gz 2_trimmed_reads/sample_R1_unpaired.fastq.gz \  
  2_trimmed_reads/sample_R2_paired.fastq.gz 2_trimmed_reads/sample_R2_unpaired.fastq.gz \  
  ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \  
  SLIDINGWINDOW:4:20 \  
  MINLEN:36
```

## 3. Quality Control on Trimmed Reads 
**Inputs:** _paired.fastq.gz  
**Outputs:** post-trimming QC reports; MultiQC summary  
**Bash:**  
```
fastqc 2_trimmed_reads/*_paired.fastq.gz -o 3_qc_trimmed/  
multiqc 3_qc_trimmed/ -o 3_qc_trimmed/
```

## 4. Generate STAR Genome Index 
**Inputs:** Genome *FASTA* file; *GTF* annotation file  
**Outputs:** *STAR* genome index files  
**Bash:**  
```
STAR --runThreadN 8 \  
  --runMode genomeGenerate \  
  --genomeDir 4_star_index/ \  
  --genomeFastaFiles genome.fa \  
  --sjdbGTFfile annotation.gtf \  
  --sjdbOverhang 99
```

## 5. Read Alignment 
**Inputs:** Trimmed paired *FASTQ* files, *STAR* genome index  
**Outputs:** *BAM* files; *STAR* log files  
**Bash:**  
```
STAR --runThreadN 8 \  
  --genomeDir 4_star_index/ \  
  --readFilesIn 2_trimmed_reads/sample_R1_paired.fastq.gz 2_trimmed_reads/sample_R2_paired.fastq.gz \  
  --readFilesCommand zcat \  
  --outFileNamePrefix 5_alignment/sample_ \  
  --outSAMtype BAM SortedByCoordinate
```

## 6. Quality Control of Aligned Files  
**Inputs:** .bam files from STAR alignment  
**Outputs:** QC summary file  
**Bash:**  
```
fastqc 5_alignment/*.bam -o 6_qc_alignment/  
multiqc 6_qc_alignment/ -o 6_qc_alignment/
```

## 7. Quantification  
**Inputs:** *BAM* files, *GTF* annotation files  
**Outputs:** Count matrix (.txt)  
**Bash:**  
```
featureCounts -T 8 -p \  
  -a annotation.gtf \  
  -o 7_featureCounts/featureCounts_matrix.txt \  
  5_alignment/*.bam
```

## 8. Differential Expression Analysis  
**Inputs:** *FeatureCounts* count matrix  
**Outputs:**
**Bash:**
