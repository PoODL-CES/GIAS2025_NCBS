Welcome to the hands-on workshop on Next-Generation Sequencing (NGS) Data Analysis. This session is designed to help you analyze Illumina-generated whole-genome resequencing (WGS) data for population genetics.
With minor adjustments, the same workflow can be applied to RAD-seq or amplicon-seq if a reference genome is available.

## Workshop Objectives: 

By the end of this workshop, you will be able to:

1. Understand and process raw sequencing data (FASTQ format)
2. Perform quality control and preprocessing
3. Map reads to a reference genome
4. Identify variants (SNPs and indels)
5. Use filtered variants for basic population genetic analysis (PCA, Admixture, heterozygosity)


## Pre-requisites & Setup

- Access to a Unix/Linux shell environment. Do the following:
  - Windows: Install [MobaXterm](https://mobaxterm.mobatek.net)  
  - Mac/Linux: Use native terminal  

- Basic knowledge of Linux commands  
  - Check the file `Understand_Linux_commands` provided in the resources  

- Conda Environment Setup  
  - To make things easier during the workshop, we have pre-installed the required conda environments for you. A file named `Conda_Environments` is included in the repository. It lists all the environments used during the sessions.

## Recommended Resources

###  Learn About Illumina Sequencing

- [How Illumina Sequencing Works](https://www.youtube.com/watch?v=fCd6B5HRaZ8&t=238s)

###  Linux Shell Basics

- Beginner: [Basic Linux Commands Exercise](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_basics.sh) \
            [Basic Linux Commands Solutions](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_basics_solutions.sh)

- Advanced :[Advanced Linux Commands Exercise](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_advanced.sh) \
            - [Advanced Linux Commands Solutions](https://github.com/PoODL-CES/Genomics_learning_workshop/blob/main/Linux_advanced_solutions.sh)


### More about Linux

- [Linux Tutorial](https://ryanstutorials.net/linuxtutorial/navigation.php)


## Workflow Overview

### 1. Raw Data to VCF

| Step | Tool(s) | Output |
|------|---------|--------|
| Quality Check | `FastQC` | QC reports |
| Trimming | `Trim-galore` | Cleaned FASTQ |
| Alignment | `BWA`| SAM/BAM |
| Post-processing | `SAMtools`, `Picard` | Sorted, deduped BAM |
| Variant Calling |  `Strelka` | Raw VCF |
| Filtering (not demonstrated here) |  `vcftools` , `bcftools`| Filtered VCF |


### 2. Population Genetics

| Step | Tool(s) | Output |
|------|---------|--------|
| Convert VCF to PLINK format (not demonstrated here) | `plink`, `vcftools` | `.bed/.bim/.fam` |
| PCA (Principal Component Analysis) | `plink` | PCA plot |
| Admixture Analysis | `ADMIXTURE` | Ancestry proportions |
