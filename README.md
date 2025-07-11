# RNA-seq Preprocessing Tutorial 
This repository has scripts for an RNA-Seq data preprocessing pipeline. It starts with raw sequencing reads and includes quality control, adapter trimming, genome indexing, alignment, and read quantification, preparing the data for downstream differential expression analysis.  

## Why Conda?
Conda is an open-source package manager and environment management system. It simplifies software installation, especially in bioinformatics where many tools are Linux-based and have complex dependencies.  
**Isolated environments:** Conda allows us to create separate environments so that different projects with different tools can run without conflict.  
**Ease of installation:** Many bioinformatics tools (ex. *STAR*, *Trimmomatic*, *FastQC*, etc) are available through [Bioconda](https://bioconda.github.io/)    

### Setting Up Conda
It's recommended to install either Miniconda (which is a lightweight version) or Anaconda (full version) for managing environments.  
**Example commands:**  
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

## Pipeline Overview 

### Adapater Trimming 
[Trimmomatic](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf) to trim adapter sequences  
Any unwanted sequences that are added when sequencing such as adapter sequences or low-quality ends of reads are removed to improve alignment accuracy.  

### Quality Control of Trimmed Files 
[FastQC](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html) on *FASTA* files  
QC before and after trimming files to see if adapter trimming improved quality of reads. This helps identify any disparities between or identify low-quality samples.  

### Generate STAR Genome Index 
[STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) on downloaded reference genome  
In order to align reads, *STAR* needs a reference genome index so that it can locate what part of the genome the RNA reads came from.  

### Read Alignment 
[STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf) to align reads to reference genome  
In this step, the reads are mapped to the reference genome, which provides the genes/regions that the RNA came from.  

### Quality Control of Aligned Files 
[MultiQC](https://docs.seqera.io/multiqc) summary report  
QC once again to see how well the reads mapped to genes and how many were counted. This produces a report of all samples, making it easier to compare between multiple samples.  

### Quantification 
[FeatureCounts](https://bioconductor.org/packages/devel/bioc/vignettes/Rsubread/inst/doc/SubreadUsersGuide.pdf) to quantify gene expression  
Here, *FeatureCounts* counts how many reads mapped to each gene. The result is a matrix where each row is a gene and each column is a sample, with the data being the number of reads for each gene in each sample.   
   

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
Download the genome file and GTF annotation files (example Rattus norvegicus files from Ensembl):  
[Rattus norvegicus genome file download](https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/dna_index/Rattus_norvegicus.GRCr8.dna.toplevel.fa.gz)  
[Rattus norvegicus GTF annotation download](https://ftp.ensembl.org/pub/release-114/gtf/rattus_norvegicus/Rattus_norvegicus.GRCr8.114.gtf.gz)  

Unzip genome and GTF files (Rattus norvegicus example): 
```
gunzip Rattus_norvegicus.mRatBN7.2.dna.primary_assembly.fa.gz
gunzip Rattus_norvegicus.mRatBN7.2.111.gtf.gz
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
