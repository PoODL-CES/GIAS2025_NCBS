### Login credentials
ssh gias3@172.16.222.186
Password:gias@ncbs 

## Navigate to your respective participant directories 

cd participant_X/ (replace X with your participant number)


### conda environment-https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

## Activate the conda environment
conda activate fastqc

#fastqc: name of the package or the tool we are installing in the new environment. Used for quality control of sequencing data. (https://mugenomicscore.missouri.edu/PDF/FastQC_Manual.pdf)

## run fastqc
fastqc input_files/*.fq.gz -o output_files

## Output files generated after running fastqc = "BEN_CI16_sub_1_fastqc.zip" and "BEN_CI16_sub_1_fastqc.html".
# .html file includes visualisations and details about the quality metrics of the sequencing data.
# .gzip file includes summary of the main quality control metrics, detailed quality control data in text format, graphs and images using in .html report.
# .gzip would have much more detailed information


## Deactivate the conda environment
conda deactivate 

## Open a new tab and download the html file to your system.

scp -r gias3@172.16.222.186:~/participant_X/output_files/*.html ~/ (X=put your participant number, type password when asked: password is gias@ncbs)

# Navigate to the directory where you have downloaded the html files in your computer and open it to visualize fastqc results



#######################################################################################################################################################################

### Trimming 
## Trimming refers to cleaning your raw sequencing reads before mapping. This step removes:adapters, low quality bases, and too short reads
## We use Trim-galore for trimming (https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md)

## activate the conda environment
conda activate trim-galore


#trim_galore: executes Trim Galore tool. it has two tools; Cutadapt (for adapter trimming) and FastQC (for QC checks)
#--paired: specifies input files are paired-end reads (DNA fragment is sequenced from both ends, producing two reads per fragment).

###  trimming all file named as fq.gz use the following command

# Use * as wildcard to consider any file with the format *_1.fq.gz as an input
for file in input_files/*_1.fq.gz
do
  mate="${file/_1.fq.gz/_2.fq.gz}"
  trim_galore --paired "$file" "$mate" --output_dir output_files
done


# "$file" refers to the current forward read (read 1).
# ${file/_1.fq.gz/_2.fq.gz} dynamically generates the corresponding reverse read (read 2) by replacing _1 with _2 in the file name.

# Output files generated after running trim-galore : 
#BEN_CI16_sub_1.fq.gz_trimming_report.txt (summary report generated during the trimming process)
#BEN_CI16_sub_1_val_1.fq.gz (val_1: refers to a validated or quality-checked file. In tools like Trim Galore, files are renamed with val_1 post trimming process)


## Deactivate the conda environment
conda deactivate 

###############################################################################################################################################################

### Mapping and sorting
# Mapping-the process of aligning sequencing reads to a reference genome to determine their original genomic location.
# Map trimmed reads using bwa mem.
# BWA (Burrows-Wheeler Aligner) — it's a fast and memory-efficient tool used to align short DNA sequencing reads to a reference genome.

## Activate the conda environment
conda activate bwa
(https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bwa/running_bwa_commands/)

#bwa mem: runs the "mem" algorithm of BWA and give SAM (Sequence Alignment/Map format) as a output file. It is optimum for 70bp-1Mbp reads, and commonly used for Illumina short-read data.
#GCA_021130815.1_PanTigT.MC.v3_genomic.fna: reference genome in FASTA format. The bwa index file should also be present.
#BEN_NW13_sub_1_val_1.fq.gz BEN_NW13_sub_2_val_2.fq.gz: These are the paired-end FASTQ files.

# Mapping all reads to reference genome in single step

for file1 in output_files/*_sub_1_val_1.fq.gz; do
    file2=${file1/_sub_1_val_1.fq.gz/_sub_2_val_2.fq.gz}
    sample_name=$(basename "$file1" _sub_1_val_1.fq.gz)

    bwa mem input_files/GCA_021130815.1_PanTigT.MC.v3_genomic.fna "$file1" "$file2" > "output_files/${sample_name}_aligned_reads.sam"
done

#for file1 in *_sub_1_val_1.fq.gz; do: looks through all the read 1(forward) FASTQ files.
#file2=${file1/_sub_1_val_1.fq.gz/_sub_2_val_2.fq.gz}: constructs the corresponding read 2 (reverse) file name by replacing _sub_1_val_1.fq.gz with _sub_2_val_2.fq.gz in file1.
#sample_name=$(basename "$file1" _sub_1_val_1.fq.gz): removes the _sub_1_val_1.fq.gz part from the filename, leaving just the sample ID (e.g., BEN_NW13).
# bwa mem GCA_021130815.1_PanTigT.MC.v3_genomic.fna "$file1" "$file2" > "${sample_name}_aligned_reads.sam": runs the actual alignment
#> BEN_NW_13_aligned_reads.sam: SAM file alignment output. It contains detailed alignment information for each read.


## Deactivate the conda environment
conda deactivate

### Convert sam to bam by samtools  (https://www.htslib.org/doc/samtools.html) and sort the read. 
# BAM-Binary Alignment/Map. Converting SAM files to BAM is an essential step for downstream analysis because it is compressed,can be indexed and more efficient for downstream analysis.  

conda activate samtools

for file in output_files/*.sam; do 
    filename=$(basename "$file" .sam)
    samtools view -S -b "$file" | samtools sort -o "output_files/${filename}_sorted.bam"
done

#samtools view: samtools view converts SAM file to BAM format.
#-S -b "$file": -s tells samtools that the input is in SAM format. -b tells that the output should be in BAM format. "$file" is the name of the current .sam file in the loop.
#samtools sort: sorts the BAM file by genomic coordinates.
#"${file%.sam}: strips off the .sam extension
#then _sorted.bam is added.
## sort: This samtools subcommand is used to sort the alignment data based on genomic coordinates.
#Sorting necessary for BAM file indexing, downstream analysis, visualization, duplicate marking.

## Deactivate the conda environment
conda deactivate 

#################################################################################################################################################################

### Mark Duplicate
# It's the process of finding and labeling copies of the same DNA  which arise during PCR amplification.
# Tools-GATK4 ( https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)

conda activate gatk4


for file in output_files/*_sorted.bam; do
    base=$(basename "$file" _sorted.bam)
    gatk MarkDuplicates \
        -I "$file" \
        -O "output_files/${base}_deduplicated.bam" \
        -M "output_files/${base}_duplication_metrics.txt" \
        --REMOVE_DUPLICATES true
done

#files would be viewed as BEN_NW_13_sorted_reads.bam
#base=${file%_sorted.bam}: ${file%_sorted.bam} removes _sorted.bam from the filename.
#-I "$file": takes the input *_sorted.bam file
#-O "${base}_deduplicated.bam" \: outputs the deduplicated BAM (with duplicates removed).
#-M "${base}_duplication_metrics.txt": writes duplication statistics to this file.
#--REMOVE_DUPLICATES true: removes the duplicate reads instead of just marking them.

## Deactivate the conda environment
conda deactivate 

######################################################################################################################################################################

### Indexing 

# indexing allows quick access to specific genomic regions and improve performance of downstream analysis tools

## Indexing the deduplicated files all at once

conda activate samtools

for bam in output_files/*_deduplicated.bam; do samtools index "$bam"; 
done

conda deactivate

## for estimating sequencing statistics like coverage per chromosome/scaffold
## Tool-qualimap (http://qualimap.conesalab.org/)

conda activate qualimap

#For generating statistics for all files 

qualimap bamqc -bam output_files/*_deduplicated.bam -outdir output_files/qualimap_results -outformat HTML

#-bam : to input bam file
#-outdir : Directory for results
#-outformat HTML :  Output in HTML 

## Navigate to the results directory
cd output_files/qualimap_results

## Open the genome_results.txt file and look for Chromosome-wise coverage and other statistics
less genome_results.txt

## Navigate back to previous directory


conda deactivate 
################################################################################################################################################
### Calling VCF (variant call format) files
# We call a VCF to identify and recoord genetic variants (like SNPs and Indels) present in a sample by comparing aligned sequencing reads to a reference genome.
# Tools-strelka (https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md)

conda activate strelka

#Strelka is run in a 2 step procedure

# Step:1 - Configuration - to specify the input data and any options pertaining to the variant calling methods themselves
# Step:2 - Workflow Execution - to specify parameters pertaining to how strelka is executed.

# For joint variant calling using all samples

    configureStrelkaGermlineWorkflow.py \
        --bam  output_files/BEN_NW13_aligned_reads_deduplicated.bam \
        --bam output_files/LGS1_aligned_reads_deduplicated.bam \
        --referenceFasta input_files/GCA_021130815.1_PanTigT.MC.v3_genomic.fna \
        --runDir output_files/strelka_output_folder                                                     # run configuration step

    output_files/strelka_output_folder/runWorkflow.py -m local -j 8                                      # workflow execution step

    #-m local: Specifies that the workflow will be executed on the local machine.
    #-j 8: Defines the number of parallel threads (CPUs) to use. 
conda deactivate


##################################################################################################################################################
