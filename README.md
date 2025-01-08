# ST4gCNV

# The ST4gCNV pipeline to use Next-Generation Sequencing (short reads) data (>=10x) to assess gene copy number varations between (homoploid) populations and closely related species.

# Step 1. Megahit for rapid assembly of NGS

##conda install -c bioconda megahit

megahit -1 input_1.fastq.gz -2 input_2.fastq.gz -o input -t 10  ###or use your customized parameters




# Step 2. Simulation of 10x reads based on genome assemblies
 
## Generate the adjacent 'windows' of megahit assembly 

samtools faidx input.final.megahit.contigs.fa

cut -f1,2 input.final.megahit.contigs.fa.fai > input.final.megahit.contigs.fa.sizes

bedtools makewindows -g  input.final.megahit.contigs.fa.sizes -w 200 -s 200 > input.final.megahit.contigs.windows200.bed

## Process the .bed file to merge the last window if it's shorter than 200bp
awk 'BEGIN {OFS="\t"} {contig=$1; start=$2; end=$3; if (contig != prev_contig) {if (NR > 1) {if (prev_end - prev_start < 200) {print prev_contig, prev_start - 200, prev_end} else {print prev_contig, prev_start, prev_end}} prev_contig=contig; prev_start=start; prev_end=end} else {if (end - start < 200) {prev_end=end} else {if (prev_start != "") {print prev_contig, prev_start, prev_end} prev_start=start; prev_end=end}}} END {print prev_contig, prev_start, prev_end}' input.final.megahit.contigs.windows200.bed > input.final.megahit.contigs.windows200_merged.bed

bedtools getfasta -fi input.final.megahit.contigs.fa -bed iput.final.megahit.contigs.windows200_merged.bed -fo input.final.megahit.contigs.windows200_merged.fasta

## create a new FASTA file where each sequence from the original is repeated 10 times##
awk '/^>/ {header=$0; next} {seq=$0; for(i=0;i<10;i++) {print header "_" i; print seq}}' input.final.megahit.contigs.windows200_merged.fasta> input.final.megahit.contigs.windows200_merged10x.fasta




# Step 3. Mapping & read coverage depth estimation on genes

bowtie2-build Reference_genome.fa Reference_genome

bowtie2 --local -x Reference_genome -f input.final.megahit.contigs.windows200_merged10x.fasta -S input.final.megahit.sam -p 20

samtools view -b input.final.megahit.sam -@ 10 -o input.final.megahit.bam

samtools sort -@ 10 -o input.final.megahit.sort.bam input.final.megahit.bam 

samtools index input.final.megahit.sort.bam



# Step 4. Gene CNV estimation based on simulated read coverages
## Simulated read depth on gene CDS###
grep "CDS" Reference_genome.gff3 > Reference_genome.cds.gff

sed 's/\S\+Parent=//g' Reference_genome.cds.gff

##coverage per base###
bedtools coverage -a Reference_genome.cds.gff  -b input.final.megahit.sort.bam -d > input.final.megahit.sort.perbase.depth

##Aggregate and Calculate Average Coverage on CDS##
awk '{ sum[$9] += $11; count[$9]++ } END { for (id in sum) print id "\t" sum[id] / count[id] }' input.final.megahit.sort.perbase.depth > input.final.megahit.sort.perCDS.depth


# Step 5. BH-corrected t-tests for detecting gene CNV differentiation between populations in R
##supposed you have n samples, merge input_1.final.megahit.sort.perCDS.depth ... input_n.final.megahit.sort.perCDS.depth into a single tab-delimit table, for example "input_01_n.perCDS.depth.txt" for downstream t-tests with Benjamini-Hochberg (BH) correction between populations###

##################bulk t-tests with BH in R platform###############

# Load the necessary library
library(dplyr)

# Set working directory (adjust this to your actual directory)
setwd("/YOUR_workDirectory/")
# Load data
gene_data <- read.table("input_01_n.perCDS.depth.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) ###simulated 10x Read depth data matrix
group1_ids <- readLines("z_NN.txt") ###sample names in population 1
group2_ids <- readLines("z_NL.txt") ###sample names in population 2


# Subset data for each group
group1_data <- gene_data[, group1_ids]
group2_data <- gene_data[, group2_ids]

# Initialize a dataframe to store the results
results <- data.frame(Gene=rownames(gene_data), T_Score=numeric(nrow(gene_data)), Mean_Group1=numeric(nrow(gene_data)), Mean_Group2=numeric(nrow(gene_data)), P_Value=numeric(nrow(gene_data)), stringsAsFactors=FALSE)


for(i in 1:nrow(gene_data)) {
  group1_values <- as.numeric(group1_data[i,])
  group2_values <- as.numeric(group2_data[i,])
  
  # Initialize a placeholder for the t-test result
  t_test_result <- list(statistic = NA, p.value = NA)

  # Attempt the t-test, using tryCatch to handle errors
  test_outcome <- tryCatch({
    t.test(group1_values, group2_values)
  }, error = function(e) {
    # If an error occurs, return NA for both statistic and p.value in the list
    return(list(statistic = NA, p.value = NA))
  })
  
  # Safely assign the t-test results to the results data frame
  results$T_Score[i] <- test_outcome$statistic
  results$P_Value[i] <- test_outcome$p.value
  
  # Calculate means for both groups
  results$Mean_Group1[i] <- mean(group1_values, na.rm = TRUE)
  results$Mean_Group2[i] <- mean(group2_values, na.rm = TRUE)
}



# Adjust p-values for multiple testing using Benjamini-Hochberg method
results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")

# Determine significance and which group is higher
results$Significance <- ifelse(results$Adjusted_P_Value < 0.05, "Yes", "No")
results$Group_Higher <- ifelse(results$Mean_Group1 > results$Mean_Group2, "Group1", "Group2")

# Select required columns to include P_Value and Adjusted_P_Value
final_results <- results[, c("Gene", "T_Score", "Mean_Group1", "Mean_Group2", "P_Value", "Adjusted_P_Value", "Significance", "Group_Higher")]


# Write the final results to a file
write.table(final_results, "input_01_n.perCDS.depth.txt.ttest.txt", sep="\t", row.names=FALSE, quote=FALSE)

##################bulk t-tests with BH in R platform END###############









