# RNA-seq Preprocessing Tutorial 
This repository has scripts for an RNA-Seq data preprocessing pipeline. It starts with raw sequencing reads and includes quality control, adapter trimming, genome indexing, alignment, and read quantification, preparing the data for downstream differential expression analysis.  

In other words, think of this pipeline as preparing your data for analysis. Each step ensures your data is high-quality, comparable, and ready for interpretation.  


## Working in the Command Line  
Many of the initial steps in an RNA-Seq pipeline are performed in a **Linux shell (bash)**. Bash allows users to execute tools, move files, and run automated workflows through scripts.  

Each script in this repo:  
- Starts with a *shebang* line (*#!/bin/bash*) so it knows to run in *bash*
- Includes *SLURM directives* (resource requests for the cluster)
- Has the Conda environment setup for the tool it uses


## Why Conda?
Conda is an open-source package manager and environment management system. It simplifies software installation, especially in bioinformatics where many tools are Linux-based and have complex dependencies.

**Isolated environments:** Conda allows us to create separate environments so that different projects with different tools can run without conflict.  
**Ease of installation:** Many bioinformatics tools (ex. *STAR*, *Trimmomatic*, *FastQC*, etc) are available through [Bioconda](https://bioconda.github.io/) , which is a community-maintained collections of bioinformatics packages.    

### Setting Up Conda
It's recommended to install either **Miniconda** (which is a lightweight version) or **Anaconda** (full version) for managing environments.  
For bioinformatics workflows (such as RNA-Seq), **Miniconda** is usually preferred since it allows users to install only the tools necessary and keeps the environment minimal.  

On the VCU Athena cluster, Anaconda and Miniconda are already preinstalled. Load either of them into your environment with: 
```
module load anaconda3/2023.09-0
```
or 
```
module load miniconda3/py39_23.9.0
```

### Basic Conda Workflow:  
Create new environment: 
```
conda create -n envName
```
Activate environment: 
```
conda activate envName
```
Install package from bioconda (example uses fastqc): 
```
conda install bioconda::fastqc
```
**It is recommended to create a new environment for each step of this pipeline. This tutorial has already included this process at each step.**


## Input Data  
RNA-Seq starts with compressed *FASTQ* files (.fastq.gz) from the sequencing facility. Each file contains the raw reads (short fragments of cDNA) plus quality scores for each base.  

For paired-end sequencing, each sample has two files:
- R1 = forward reads (Sample1_R1.fastq.gz)  
- R2 = reverse reads (Sample1_R2.fastq.gz)  

R1 and R2 together represent both ends of each fragment. These paired *FASTQ* files are the starting point of the preprocessing pipeline.


## Pipeline Overview 

### Adapater Trimming 
[Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) to trim adapter sequences  
Any unwanted sequences that are added when sequencing such as adapter sequences or low-quality ends of reads are removed to improve alignment accuracy and overall quality of the reads. In order to perform this step accurately, users need to know the specific adapter sequences that were used during library preparation. These are usually provided by the sequencing facility.  

### Merging Paired-End Reads (Optional)  
[BBMerge](https://anaconda.org/agbiome/bbtools) to merge paired-end reads  
Paired-end reads often overlap and BBMerge identifies these regions to combine them into a single, higher-quality consensus read. Reads that cannot be merged remain as separate forward (R1) and reverse (R2) files. This step is not always necessary but can simplify downstream analysis when insert sizes are short.

### Quality Control of Trimmed Files 
[FastQC](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html) on *FASTA* files  
FastQC is used to assess the quality of sequencing reads both **before** and **after** trimming. Running FastQC on raw files provides a baseline quality report, while repeating it after trimming confirms that adapters and low-quality bases were successfully removed. These reports can also reveal problematic samples early, saving time in downstream analysis.  

### Generate STAR Genome Index 
[STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) on downloaded reference genome  
Before aligning reads, *STAR* requires a reference genome index that tells it where genes and exons are located. To build this index, you’ll need:  
- A reference genome in *FASTA* format (.fa)
- A gene annotation file in *GTF* or *GFF3* format (.gtf, .gff3)

Both files must come from the same genome build (Ensembl rn6). Be sure to decompress .gz files before use. The index is built once and can be reused for multiple samples.  

### Read Alignment 
[STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) to align reads to reference genome  
Reads are mapped to the reference genome in a process called alignment (or mapping), which determines where each RNA-Seq read originated. *STAR* is optimized for RNA-Seq, handling spliced reads that span exon–exon junctions. The output includes *BAM* files (aligned reads) and log files with alignment statistics. Accurate alignment is crucial for reliable quantification.  

### Quality Control of Aligned Files 
[MultiQC](https://docs.seqera.io/multiqc) summary report  
After alignment, QC is performed again to evaluate mapping quality across all samples. *MultiQC* compiles the results from multiple tools (FastQC, STAR logs) into one combined report. This summary makes it easier to compare multiple samples side by side and quickly identify inconsistencies.  

### Quantification 
[FeatureCounts](https://bioconductor.org/packages/devel/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) to quantify gene expression  
After alignment, *featureCounts* counts how many reads overlap with annotated genes. The output is a **count matrix** where each row corresponds to a gene and each column corresponds to a sample. This table is the foundation for downstream analyses (differential expression using DESeq2 or edgeR). Quantification can be performed entirely on the command line or within R.   


## Running Jobs on the Cluster with Slurm    
When working on a **high-performance computing (HPC) cluster**, we don’t run heavy analyses directly on the login node. The login node is shared by all users and is only meant for light tasks, like preparing scripts, moving files, submitting jobs.

All compute-intensive steps (like trimming, aligning, and quantifying RNA-Seq data) are run on **compute nodes**. To access these, you must submit jobs through the cluster’s scheduler.

### What is SLURM?  
Our cluster uses **SLURM (Simple Linux Utility for Resource Management)**, which schedules and manages jobs on compute nodes. **SLURM** decides when and where your job runs based on available resources.  

### Job Scripts  
To run analyses on compute nodes, you create a job script:
- Starts with a *shebang* line (*#!/bin/bash*) to ensure *Bash* is used. Without it, the system may default to another shell, causing errors
- Contains *SLURM directives* (*#SBATCH*) that request resources (CPUs, memory, time)
- Includes the commands to run your tool or pipeline step

Each step of this RNA-Seq pipeline has its own script in the *scripts/* folder, already formatted with a *shebang* line and *SLURM directives*.

### Submitting Jobs  
**1. Make the script executable:**  
```
chmod +x script.sh
```
The filename will often turn green in your terminal, which means it’s executable.

**2. Submit the script to Slurm:**  
```
sbatch script.sh
```

**3. Check job status:**  
```
squeue -u username
```
or  
```
squeue #jobnumber
```

## 1. Setting Up 
First, be sure to set your current directory to wherever your scripts are located (depending on which script you are running): 
```
cd *path to scripts* 
```
Next, if you are working in a cluster, make sure to load the anaconda module (adjust accordingly to which version you are using):

```
module load anaconda3/2023.09-0
```
Make sure to do these steps everytime you restart the shell so that everything is loaded.  


## 2. Adapter Trimming 
**Inputs:** _R1.fastq.gz & _R2.fastq.gz  
**Outputs:** Paired and unpaired trimmed *FASTQ* files  
### 2a. Quality Check with FastQC
**Bash:**  
Create environment to run FastQC (fqcEnv): 
```
conda create -n fqcEnv
```
Activate environment: 
```
conda activate fqcEnv
```
Install FastQC in environment from bioconda: 
```
conda install bioconda::fastqc
```
Activate script for FastQC (FastQC.sh): 
```
chmod +x FastQC.sh
```
Submit job for FastQC script: 
```
sbatch FastQC.sh
```
### 2b. Trim Reads  
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
### 2c. Quality Check on trimmed reads
Rerun FastQC using steps from **2a** on trimmed files to see whether quality of reads improved after trimming adapter sequences.  


## 3. Generate STAR Genome Index 
**Inputs:** Genome *FASTA* file; *GTF* annotation file  
**Outputs:** *STAR* genome index files  
**Bash:**  
Download the genome file and GTF annotation files (example Mus musculus files from Ensembl):  
[Mus musculus genome file download](https://ftp.ensembl.org/pub/release-115/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.toplevel.fa.gz)  
[Mus musculus GTF annotation download](https://ftp.ensembl.org/pub/release-115/gtf/mus_musculus/Mus_musculus.GRCm39.115.gtf.gz)  

Unzip genome and GTF files (Mus musculus example): 
```
gunzip Mus_musculus.GRCm39.dna.toplevel.fa.gz
gunzip Mus_musculus.GRCm39.115.gtf.gz
```
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
Make sure to have the STARallsamples script saved as well because that is what we use to loop through all given samples!  
**Bash:**  
Activate STAR environment (if deactivated): 
```
conda activate starEnv
```
Activate script for STARaligner (STARaligner.sh): 
```
chmod +x STARaligner.sh
```
Submit job for STARaligner script: 
```
sbatch STARaligner.sh
```


## 5. Quality Control of Aligned Files  
**Inputs:** .bam files from STAR alignment  
**Outputs:** MultiQC summary file  
**Bash:**  
Create environment to run MultiQC (mqcEnv): 
```
conda create -n mqcEnv
```
Activate environment: 
```
conda activate mqcEnv
```
Install MultiQC in environment from bioconda: 
```
conda install bioconda::multiqc
```
Activate script for MultiQC (MultiQC.sh): 
```
chmod +x MultiQC.sh
```
Submit job for MultiQC script: 
```
sbatch MultiQC.sh
```


## 6. Quantification  
**Inputs:** *BAM* files, *GTF* annotation files  
**Outputs:** Count matrix (.txt)  
**Bash:**  
Create environment to run featureCounts (subreadEnv): 
```
conda create -n subreadEnv
```
Activate environment: 
```
conda activate subreadEnv
```
Install Subread in environment from bioconda: 
```
conda install bioconda::subread
```
Activate script for featureCounts (featureCounts.sh): 
```
chmod +x featureCounts.sh
```
Submit job for featureCounts script: 
```
sbatch featureCounts.sh
```
