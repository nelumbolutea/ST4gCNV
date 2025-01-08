##################bulk t-tests with BH in R platform###############

## Load the necessary library
library(dplyr)

## Set working directory (adjust this to your actual directory)
setwd("/YOUR_workDirectory/")
## Load data
gene_data <- read.table("input_01_n.perCDS.depth.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE) ###simulated 10x Read depth data matrix
group1_ids <- readLines("z_NN.txt") ###sample names in population 1
group2_ids <- readLines("z_NL.txt") ###sample names in population 2

## Subset data for each group
group1_data <- gene_data[, group1_ids]
group2_data <- gene_data[, group2_ids]

### Initialize a dataframe to store the results
results <- data.frame(Gene=rownames(gene_data), T_Score=numeric(nrow(gene_data)), Mean_Group1=numeric(nrow(gene_data)), Mean_Group2=numeric(nrow(gene_data)), P_Value=numeric(nrow(gene_data)), stringsAsFactors=FALSE)

for(i in 1:nrow(gene_data)) {
  group1_values <- as.numeric(group1_data[i,])
  group2_values <- as.numeric(group2_data[i,])
  
  ### Initialize a placeholder for the t-test result
  t_test_result <- list(statistic = NA, p.value = NA)

  ### Attempt the t-test, using tryCatch to handle errors
  test_outcome <- tryCatch({
    t.test(group1_values, group2_values)
  }, error = function(e) {
    # If an error occurs, return NA for both statistic and p.value in the list
    return(list(statistic = NA, p.value = NA))
  })
  
  ### Safely assign the t-test results to the results data frame
  results$T_Score[i] <- test_outcome$statistic
  results$P_Value[i] <- test_outcome$p.value
  
  ### Calculate means for both groups
  results$Mean_Group1[i] <- mean(group1_values, na.rm = TRUE)
  results$Mean_Group2[i] <- mean(group2_values, na.rm = TRUE)
}



## Adjust p-values for multiple testing using Benjamini-Hochberg method
results$Adjusted_P_Value <- p.adjust(results$P_Value, method = "BH")

## Determine significance and which group is higher
results$Significance <- ifelse(results$Adjusted_P_Value < 0.05, "Yes", "No")
results$Group_Higher <- ifelse(results$Mean_Group1 > results$Mean_Group2, "Group1", "Group2")

## Select required columns to include P_Value and Adjusted_P_Value
final_results <- results[, c("Gene", "T_Score", "Mean_Group1", "Mean_Group2", "P_Value", "Adjusted_P_Value", "Significance", "Group_Higher")]


## Write the final results to a file
write.table(final_results, "input_01_n.perCDS.depth.txt.ttest.txt", sep="\t", row.names=FALSE, quote=FALSE)

##################bulk t-tests with BH in R platform END###############