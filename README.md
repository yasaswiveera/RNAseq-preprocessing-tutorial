# RNA-seq Tutorial 
This repository contains scripts for a full RNA-Seq analysis pipeline. Starting with raw sequencing reads, it undergoes quality control, read trimming, alignment, quantification, and statistical analysis to obtain differential expression results. 

## Pipeline Overview 

### Adapater Trimming 
**Trimmomatic** to trim adapter sequences  
Any unwanted sequences that are added when sequencing such as adapter sequences or low-quality ends of reads are removed to improve alignment accuracy.  

### Generate STAR Genome Index 
**STAR** on downloaded reference genome  
In order to align reads, *STAR* needs a reference genome index so that it can locate what part of the genome the RNA reads came from.  

### Read Alignment 
**STAR** to align reads to reference genome  
In this step, the reads are mapped to the reference genome, which provides the genes/regions that the RNA came from.  

### Quantification 
**FeatureCounts** to quantify gene expression  
Here, *FeatureCounts* counts how many reads mapped to each gene. The result is a matrix where each row is a gene and each column is a sample, with the data being the number of reads for each gene in each sample.   

### Quality Control of Aligned Files 
**FastQC** on *BAM* files  
**MultiQC** summary report  
QC once again to see how well the reads mapped to genes and how many were counted. This further helps identify any disparities between or identify low-quality samples.   

### Differential Expression Analysis 
**EdgeR** for statistical analysis; includes data normalization, exploratory analyses, and visualization techniques  
Here, statistical modeling is used to identify differentially expressed genes between experimental groups. Further, *EdgeR* can help visualize any trends between samples and genes.    

## 1. Setting Up 
First, be sure to set your current directory to wherever your scripts are located (depending on which script you are running): 
```
cd *path to scripts* 
```
Next, if you are working in a cluster, make sure to load the anaconda module (adjust accordingly to which version you are using):

```
module load anaconda3/2023.09-0
```

## 2. Adapter Trimming 
**Inputs:** _R1.fastq.gz & _R2.fastq.gz  
**Outputs:** Paired and unpaired trimmed *FASTQ* files  
**Bash:**  
Create environment to run Trimmomatic (trimEnv) and make sure to update this name in the script accordingly: 
```
conda create -n trimEnv
```
Activate environment: 
```
conda activate trimEnv
```
Install Trimmomatic in environment from bioconda: 
```
conda install bioconda::trimmomatic
```
Activate script for Trimmomatic (trimmomatic.sh): 
```
chmod +x trimmomatic.sh
```
Submit job for Trimmomatic script: 
```
sbatch trimmomatic.sh
```

## 3. Generate STAR Genome Index 
**Inputs:** Genome *FASTA* file; *GTF* annotation file  
**Outputs:** *STAR* genome index files  
**Bash:**  
Download the genome file and GTF annotation files (example Rattus norvegicus files from Ensembl):  
[Rattus norvegicus genome file download](https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/dna_index/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz)
[Rattus norvegicus GTF annotation download](https://ftp.ensembl.org/pub/release-114/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.114.gtf.gz)  

Create environment to run STAR (starEnv): 
```
conda create -n starEnv
```
Activate environment: 
```
conda activate starEnv
```
Install STAR in environment from bioconda: 
```
conda install bioconda::star
```
Activate script for STARindex (STARindex.sh): 
```
chmod +x STARindex.sh
```
Submit job for STARindex script: 
```
sbatch STARindex.sh
```

## 4. Read Alignment 
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

## 5. Quality Control of Aligned Files  
**Inputs:** .bam files from STAR alignment  
**Outputs:** QC summary file  
**Bash:**  
```
fastqc 5_alignment/*.bam -o 6_qc_alignment/  
multiqc 6_qc_alignment/ -o 6_qc_alignment/
```

## 6. Quantification  
**Inputs:** *BAM* files, *GTF* annotation files  
**Outputs:** Count matrix (.txt)  
**Bash:**  
```
featureCounts -T 8 -p \  
  -a annotation.gtf \  
  -o 7_featureCounts/featureCounts_matrix.txt \  
  5_alignment/*.bam
```

## 7. Differential Expression Analysis  
**Inputs:** *FeatureCounts* count matrix  
**Outputs:**
**Bash:**
