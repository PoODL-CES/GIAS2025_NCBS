#PCA

#firstly, .bed, .bim and .fam files need to be generated. from a vcf file. 
#.bed: binary file that stores genotype calls for every individual at very SNP.
#.bim: tab-delimited text file with information about each SNP.
#.fam: tab-delimited text file; contains information about each individual/sample.

 conda activate plink
 plink --bfile output_file \
  --pca 10 \
  --allow-extra-chr \
  --out /home/gias3/admin/output_files/output_file_pca

#--bfile output_file: loads the plink binary dataset previously created using --make-bed
#--pca 10: calculates the top 10 principal components
#--allow-extra-chr: allows non-standard chromosome names
#out output_file_pca: sets the prefix for output files (example: output_file_pca.eigenval)
#eigenvalues and eigenvectors would be created

conda activate r_env
R
install.packages("ggplot2")
library(ggplot2)      
eigenvec_data <- read.table("output_file_pca.eigenvec", header=FALSE)
colnames(eigenvec_data) <- c("FID", "IID", paste("PC", 1:10, sep=""))
head(eigenvec_data)
pdf(file = paste0(output_dir, "pca_plot.pdf"))
ggplot(eigenvec_data, aes(x=PC1, y=PC2)) +
  geom_point() +
  labs(x="Principal Component 1", y="Principal Component 2", title="PCA Plot: PC1 vs PC2") +
  theme_minimal()
  dev.off()
ggsave("/home/gias3/admin/output_files/pca_plot.png")
q()


scp username@IP_address:~/"path to the file on the remote cluster"/Rplots.pdf .
scp username@IP_address:~/"path to the file on the remote cluster"/pca_plot.png .

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

INPUT_DIR="/home/gias3/admin/input_files"
OUTPUT_DIR="/home/gias3/admin/output_files"
INPUT_PREFIX="machali_Aligned_rangeWideMerge_strelka_update2_BENGAL_mac3_passOnly_biallelicOnly_noIndels_minMAF0Pt05_chr_E2_minDP3_minQ30_minGQ30_hwe_0.05_noIndels_missing_mm0.6_meandepth95percentile_noZSB"

mkdir -p "$output_files"

for K in {2..4}
do 
  admixture "$INPUT_DIR/$INPUT_PREFIX.bed" $K --out "$OUTPUT_DIR/${INPUT_PREFIX}_K${K}"
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
