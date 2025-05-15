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
