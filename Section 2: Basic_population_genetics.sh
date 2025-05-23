### Login credentials
ssh gias3@172.16.222.186
Password:gias@ncbs 

#PCA

#firstly, .bed, .bim and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.

conda activate plink

plink --bfile input_files/machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB \
  --pca 5 \
  --allow-extra-chr \
  --out machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB

conda deactivate

#--pca 5: calculates the top 5 principal components
#--allow-extra-chr: allows non-standard chromosome names
#out output_file_pca: sets the prefix for output files (example: output_file_pca.eigenval)
#eigenvalues and eigenvectors would be created

#to plot the PCA we will use R (https://www.r-project.org)
conda activate r-ggplot2 

R
library(ggplot2)
# Read in the PCA data
eigenvec_data <- read.table("machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.eigenvec", header=FALSE)

# Add column names (FID, IID, Geographic Region, PCs)
colnames(eigenvec_data) <- c("FID", "IID", paste0("PC", 1:5))

# Plot
ggplot(eigenvec_data, aes(x=PC1, y=PC2, color=IID)) +
  geom_point(size=3, alpha=0.8) +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_color_brewer(palette = "Set2") # You can change to "Dark2", "Paired", etc.

# Save plots
ggsave("/home/gias3/admin/output_files/pca_plot.pdf")
ggsave("/home/gias3/admin/output_files/pca_plot.png")

q()
#Save workspace image? [y/n/c]: y
pwd (#print working directory)
exit (#exit from the GIAS cluster)

scp username@IP_address:~/"path to the file on the remote cluster"/Rplots.pdf .
scp username@IP_address:~/"path to the file on the remote cluster"/pca_plot.png .
#example: scp gias3@172.16.222.186:/home/gias3/participant_3/output_files/Rplots.pdf .

#R: Launches R
#install.packages("ggplot2") library(ggplot2): installs and loads ggplot2 for plotting
#eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE): reads the .eigenvec file into the dataframe
#colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep="")): assigns column names; first 2 columns are FID and IID while next 10 columns are PC1 to PC10
#head(eigenvec_data): displays the first few rows for the purpose of confirmation
#ggplot(eigenvec_data, aes(x=PC1, y=PC2)) + geom_point() + labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") + theme_minimal()
       # Actually builds the graph
#ggsave("pca_plot.png"): saves the plot in .png format
#q(): Exits the R console

#scp: secure copy protocol; copies files betweeen a local and a remote computer.
# .: destination on the local machine. copies to the current directory.
#pca_plot.png and Rplots.pdf would be saved in your local home directory


#Admixture

conda activate admixture

INPUT_DIR="/home/gias3/participant_X/input_files"
OUTPUT_DIR="/home/gias3/participant_X/output_files"
INPUT_PREFIX="machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB"

for K in {2..4}
do
  admixture "$INPUT_DIR/$INPUT_PREFIX.bed" $K
  mv "${INPUT_PREFIX}.${K}.Q" "$OUTPUT_DIR/${INPUT_PREFIX}_K${K}.Q"
  mv "${INPUT_PREFIX}.${K}.P" "$OUTPUT_DIR/${INPUT_PREFIX}_K${K}.P"
done

#for K in {2..4}: loops through the values 2, 3, and 4
#INPUT_DIR: specifies the location of input directory
#OUTPUT_DIR: specifies the location of output directory
#INPUT_PREFIX: base name of the PLINK file (without the .bed/.bim/.fam extension).
#for K in {2..4}: runs the simulation for 2 to 4 ancestral populations
#${K}: number of ancestral populations we are trying to estimate (2 to 4 in this case)

#output would be as follows:

#for K=2;

#input_file_K2.Q: Q matrix (individual ancestry populations)
#input_file_K2.P: P matrix (ancestral allele frequencies)
#input_file_K2.log: Log file with convergence and likelihood info

#same pattern for K=3 and K=4

#To visualise the ADMIXTURE results, we plot it in R

conda activate r_env

R #(It will open R prompt)

# Load libraries 

library(ggplot2) 

library(reshape2) 

library(dplyr)  

# Read ADMIXTURE Q file and sample IDs from .fam 

q3 <- read.table("machali_...3.Q")  # Replace with actual full filename 

fam <- read.table("../input_files/machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.fam")

sample_ids <- fam$V1 

# Add sample IDs 

q3$ID <- sample_ids 


# OPTIONAL: Add groups/populations (e.g., from file or by pattern) 

# If you have a file with sample ID and population: 

# groups <- read.table("group_info.txt", header = TRUE) 

# q3 <- merge(q3, groups, by.x = "ID", by.y = "SampleID") 

# OR, extract group info from sample ID 

q3$Group <- gsub(".*_(.*)", "\\1", q3$ID)  # example: get last part after underscore 

# Reshape data for ggplot 

q3_long <- melt(q3, id.vars = c("ID", "Group")) 

# Sort individuals by group 

q3_long <- q3_long %>% 

  arrange(Group, ID) %>% 

  mutate(ID = factor(ID, levels = unique(ID))) 

# Plot 

p <- ggplot(q3_long, aes(x = ID, y = value, fill = variable)) + 

  geom_bar(stat = "identity", width = 1) + 

  facet_grid(~Group, scales = "free_x", space = "free_x") + 

  theme_minimal() + 

  labs(x = "Individuals", y = "Ancestry Proportion", title = "ADMIXTURE Plot (K=3)") + 

  theme(axis.text.x = element_blank(), 

        axis.ticks.x = element_blank(), 

        panel.spacing = unit(0.5, "lines"), 

        strip.text.x = element_text(angle = 0, face = "bold"), 

        legend.position = "right") + 

  scale_fill_brewer(palette = "Set1") 
# Save as PNG 

ggsave("admixture_K3_grouped.png", p, width = 12, height = 6, dpi = 300) 

#For Heterozygosity
conda activate rtg-tools
rtg vcfstats input_files/machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.recode.vcf.gz > output_files/machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.vcfstats
conda deactivate

less machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB.vcfstats

#rtg vcfstats ... > /home/gias3/admin/output_files/... : Run vcfstats to generate statistics from the VCF file and save the results to the specified output directory
