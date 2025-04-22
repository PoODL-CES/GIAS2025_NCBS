Welcome to the hands-on workshop on Next-Generation Sequencing (NGS) Data Analysis. This session is designed to help you analyze Illumina-generated whole-genome resequencing (WGS) data for population genetics.
With minor adjustments, the same workflow can be applied to RAD-seq or amplicon-seq if a reference genome is available.

## Workshop Objectives: 

By the end of this workshop, you will be able to:

1. Understand and process raw sequencing data (FASTQ format)
2. Perform quality control and preprocessing
3. Map reads to a reference genome
4. Call and filter variants (SNPs and indels)
5. Use filtered VCFs for population structure analysis (PCA, Admixture)


## Pre-requisites & Setup

- Access to a Unix/Linux shell environment
  - Windows: Install [MobaXterm](https://mobaxterm.mobatek.net)
  - Mac/Linux: Use native terminal

- Basic knowledge of Linux commands
  - Check the file `Understand Linux commands` provided in the resources

## Recommended Resources

###  Learn About Sequencing

- [How Illumina Sequencing Works](https://www.youtube.com/watch?v=fCd6B5HRaZ8&t=238s)

###  Linux Shell Basics

- Beginner: `Linux_basics.sh`, `Linux_basics_solutions.sh`
- Advanced: `Linux_advanced.sh`, `Linux_advanced_solutions.sh`


## Workflow Overview

### 1. Raw Data to VCF

| Step | Tool(s) | Output |
|------|---------|--------|
| Quality Check | `FastQC` | QC reports |
| Trimming (optional) | `Trim-galore`, `Trimmomatic` | Cleaned FASTQ |
| Alignment | `BWA`, `Bowtie2` | SAM/BAM |
| Post-processing | `SAMtools`, `Picard` | Sorted, deduped BAM |
| Variant Calling | `bcftools`, `GATK` | Raw VCF |
| Filtering | `bcftools`, `vcftools` | Filtered VCF |

### 2. VCF to Population Genetics

| Step | Tool(s) | Output |
|------|---------|--------|
| Convert VCF to PLINK format | `plink`, `vcftools` | `.bed/.bim/.fam` |
| PCA (Principal Component Analysis) | `plink` | PCA plot |
| Admixture Analysis | `ADMIXTURE` | Ancestry proportions |
