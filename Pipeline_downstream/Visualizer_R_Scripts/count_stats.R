library(readxl)
library(openxlsx)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)

#  =========================> <=================================
# This function grabs the title from the file and adds them to the output
read_and_rename_gene <- function(file_path) {
  # Extract sample label from the file name
  sample_label <- gsub("(.*)-DECalls-FC_gene.txt", "\\1", basename(file_path))
  
  # Read the data file
  data <- read_tsv(file_path, col_types = cols()) %>%
    rename_at(vars(log2FC, FC), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(TPM_1, TPM_2), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(Bin),~paste(sample_label, ., sep = " "))
  return(data)
}

read_and_rename_isoform <- function(file_path) {
  # Extract sample label from the file name
  sample_label <- gsub("(.*)-DECalls-FC_isoform.txt", "\\1", basename(file_path))
  
  # Read the data file
  data <- read_tsv(file_path, col_types = cols()) %>%
    rename_at(vars(log2FC, FC), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(TPM_1, TPM_2), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(Bin),~paste(sample_label, ., sep = " "))
  return(data)
}

# Read and rename for each dataset
data1 <- read_and_rename_gene("IsoDE_raw_C42/37-38-DECalls-FC_gene.txt")
data2 <- read_and_rename_gene("IsoDE_raw_C42/37-39-DECalls-FC_gene.txt")
data3 <- read_and_rename_gene("IsoDE_raw_C42/38-39-DECalls-FC_gene.txt")
# Combine the data
combined <- cbind(data1, data2, data3)
# Save combined data to an Excel file
write.xlsx(combined, file = "Binning_stats_C42_Genes.xlsx")

# now do the same for isoforms
data1a <- read_and_rename_isoform("IsoDE_raw_C42/37-38-DECalls-FC_isoform.txt")
data2a <- read_and_rename_isoform("IsoDE_raw_C42/37-39-DECalls-FC_isoform.txt")
data3a <- read_and_rename_isoform("IsoDE_raw_C42/38-39-DECalls-FC_isoform.txt")
combined_iso <- cbind(data1a, data2a, data3a)
# Save combined data to an Excel file
write.xlsx(combined_iso, file = "Binning_stats_C42_isoforms.xlsx")

#===============================================================================
#For LnCap
data4 <- read_and_rename_gene("IsoDE_raw_LnCap/46-47-DECalls-FC_gene.txt")
data5 <- read_and_rename_gene("IsoDE_raw_LnCap/46-48-DECalls-FC_gene.txt")
data6 <- read_and_rename_gene("IsoDE_raw_LnCap/47-48-DECalls-FC_gene.txt")

combined_LnCap <- cbind(data4, data5, data6)
write.xlsx(combined_LnCap, file = "Binning_stats_LnCap_Genes.xlsx")
# now do the same for isoforms
data4a <- read_and_rename_isoform("IsoDE_raw_LnCap/46-47-DECalls-FC_isoform.txt")
data5a <- read_and_rename_isoform("IsoDE_raw_LnCap/46-48-DECalls-FC_isoform.txt")
data6a <- read_and_rename_isoform("IsoDE_raw_LnCap/47-48-DECalls-FC_isoform.txt")

combined_iso_ln <- cbind(data4a, data5a, data6a)
# Save combined data to an Excel file
write.xlsx(combined_iso_ln, file = "Binning_stats_LnCap_isoforms.xlsx")

# Read the data from the text file, specifying tab as the delimiter

# Print the structure of the data frame
str(data)


write.xlsx(combined_LnCap, file = "Binning_stats_LnCap_Genes.xlsx")
# ------- Now, read the excel and plot! -------------------
db_c <- read_excel("Binning_stats_C42_Genes.xlsx")
db_l <- read_excel("Binning_stats_LnCap_Genes.xlsx")
db_l_iso <- read_excel("Binning_stats_LnCap_isoforms.xlsx")
db_c_iso <- read_excel("Binning_stats_C42_isoforms.xlsx")

# Define a function to count satisfying and not satisfying conditions
count_conditions <- function(db, prefix, orth_name1, orth_name2) {
  # Construct column names dynamically based on the prefix
  col_tpm_1 <- paste(prefix, "TPM_1", sep = " ")
  col_tpm_2 <- paste(prefix, "TPM_2", sep = " ")
  col_bin <- paste(prefix, "Bin", sep = " ")
  
  # Filter rows satisfying the condition
  filtered_db <- db %>%
    filter((.[[col_tpm_1]] >= 0.5 | .[[col_tpm_2]] >= 0.5))
  
  # Count the number of rows satisfying the condition
  num_satisfying_condition <- nrow(filtered_db)
  num_not_satisfying_condition <- nrow(db) - num_satisfying_condition
  
  # Count how many of these have specific 'Bin' values
  num_OR_1 <- sum(filtered_db[[col_bin]] == orth_name1)
  num_OR_2 <- sum(filtered_db[[col_bin]] == orth_name2)
  
  # Print results
  cat("Total samples satisfying the condition for", prefix, ":",
      num_satisfying_condition, "\n")
  cat("Number of samples not satisfying the condition for", prefix, ":",
      num_not_satisfying_condition, "\n")
  cat("Samples with",orth_name1," in '", col_bin, "' for", prefix, ":", num_OR_1, "\n")
  cat("Samples with",orth_name2, "in '", col_bin, "' for", prefix, ":", num_OR_2, "\n")
}

#get the data for stat purposes.
count_conditions(db_c, "37-38","OR_37", "OR_38")
count_conditions(db_c, "38-39","OR_38", "OR_39")
count_conditions(db_c, "37-39", "OR_37", "OR_39")

count_conditions(db_l, "46-47","OR_46", "OR_47")
count_conditions(db_l, "46-48","OR_46", "OR_48")
count_conditions(db_l, "47-48", "OR_47", "OR_48")

# isoforms:
count_conditions(db_c_iso, "37-38","OR_37", "OR_38")
count_conditions(db_c_iso, "38-39","OR_38", "OR_39")
count_conditions(db_c_iso, "37-39", "OR_37", "OR_39")

count_conditions(db_l_iso, "46-47","OR_46", "OR_47")
count_conditions(db_l_iso, "46-48","OR_46", "OR_48")
count_conditions(db_l_iso, "47-48", "OR_47", "OR_48")

