library(readxl)
library(openxlsx)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)

#  =============================================================================
#  Author: Xiaoyue Deng                       Contact: xiaoyued.deng@gmail.com
#  =============================================================================
# This function grabs the title from the file and adds them to the output
read_and_rename_gene <- function(file_path) {
  # Extract sample label from the file name
  sample_label <- gsub("(.*)-DECalls-FC_gene.txt",
                       "\\1", basename(file_path))
  
  # Read the data file
  data <- read_tsv(file_path, col_types = cols()) %>%
    rename_at(vars(log2FC, FC), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(TPM_1, TPM_2), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(Bin),~paste(sample_label, ., sep = " "))
  return(data)
}
# essentially same function, but for isoforms
read_and_rename_isoform <- function(file_path) {
  # Extract sample label from the file name
  sample_label <- gsub("(.*)-DECalls-FC_isoform.txt",
                       "\\1", basename(file_path))
  
  # Read the data file
  data <- read_tsv(file_path, col_types = cols()) %>%
    rename_at(vars(log2FC, FC), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(TPM_1, TPM_2), ~paste(sample_label, ., sep = " ")) %>%
    rename_at(vars(Bin),~paste(sample_label, ., sep = " "))
  return(data)
}
#===============================================================================
# Read and rename for each dataset
data1 <- read_and_rename_gene("IsoDE_raw_C42/37-38-DECalls-FC_gene.txt")
data2 <- read_and_rename_gene("IsoDE_raw_C42/37-39-DECalls-FC_gene.txt")
data3 <- read_and_rename_gene("IsoDE_raw_C42/38-39-DECalls-FC_gene.txt")
# now do the same for isoforms
data1a <- read_and_rename_isoform("IsoDE_raw_C42/37-38-DECalls-FC_isoform.txt")
data2a <- read_and_rename_isoform("IsoDE_raw_C42/37-39-DECalls-FC_isoform.txt")
data3a <- read_and_rename_isoform("IsoDE_raw_C42/38-39-DECalls-FC_isoform.txt")
# --------------- Export files to a collective excel ---------------------------
combined_iso <- cbind(data1a, data2a, data3a)
combined <- cbind(data1, data2, data3)
# Save combined data to an Excel file
write.xlsx(combined, file = "Binning_stats_C42_Genes.xlsx")
# Save combined data to an Excel file
write.xlsx(combined_iso, file = "Binning_stats_C42_isoforms.xlsx")
#===============================================================================
#For LnCap
data4 <- read_and_rename_gene("IsoDE_raw_LnCap/46-47-DECalls-FC_gene.txt")
data5 <- read_and_rename_gene("IsoDE_raw_LnCap/46-48-DECalls-FC_gene.txt")
data6 <- read_and_rename_gene("IsoDE_raw_LnCap/47-48-DECalls-FC_gene.txt")
# now do the same for isoforms
data4a <- read_and_rename_isoform("IsoDE_raw_LnCap/46-47-DECalls-FC_isoform.txt")
data5a <- read_and_rename_isoform("IsoDE_raw_LnCap/46-48-DECalls-FC_isoform.txt")
data6a <- read_and_rename_isoform("IsoDE_raw_LnCap/47-48-DECalls-FC_isoform.txt")
#Combine data and export as an excel
combined_LnCap <- cbind(data4, data5, data6)
write.xlsx(combined_LnCap, file = "Binning_stats_LnCap_Genes.xlsx")
combined_iso_ln <- cbind(data4a, data5a, data6a)
write.xlsx(combined_iso_ln, file = "Binning_stats_LnCap_isoforms.xlsx")

#==============================Plotting session=================================
#Load the collective databases
db_c <- read_excel("Binning_stats_C42_Genes.xlsx")
db_l <- read_excel("Binning_stats_LnCap_Genes.xlsx")
db_l_iso <- read_excel("Binning_stats_LnCap_isoforms.xlsx")
db_c_iso <- read_excel("Binning_stats_C42_isoforms.xlsx")
# ================> Function to generate rocket plot <==========================
create_scatter_plot <- function(db, prefix, ortho_name1, ortho_name2,
                                condition1, condition2) {
  # !!! Ideally: add an option about whether to save the png or not.
  # Construct column names dynamically based on the prefix
  col_tpm_1 <- paste(prefix, "TPM_1", sep = " ")
  col_tpm_2 <- paste(prefix, "TPM_2", sep = " ")
  col_bin <- paste(prefix, "Bin", sep = " ")
  
  # Filter and select data for plotting
  filtered_db <- db %>%
    filter((
      !!rlang::sym(col_tpm_1) >= 0.5) | (!!rlang::sym(col_tpm_2) >= 0.5)) %>%
    select(all_of(c(col_tpm_1, col_tpm_2, col_bin)))
  
  # Create a scatter plot
  plot <- ggplot(filtered_db,
                 aes(x = log10(get(col_tpm_2)), y = log10(get(col_tpm_1)),
                                  color = get(col_bin))) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = setNames(c("#a8e6cf", "#ff8b94"),
                                         c(ortho_name1, ortho_name2))) +
    labs(title = paste("Scatter Plot of Log10 Transformed TPM Values (",
                       condition1, " vs ", condition2, ") , genes"),
         x = paste("Log10(", condition2, ")"),
         y = paste("Log10(", condition1, ")"),
         color = "Gene Bin") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white")) +
    theme(legend.position = "right")
  # Save the plot
  file_name <- paste("Log10_TPM_Scatter_Plot_", condition1,"vs",
                     condition2 , ".png", sep="")
  ggsave(file_name, plot, width = 8, height = 6)
  # Print the plot
  print(plot)
}
# ===========================> Usage example <==================================
#===========================> Combined Plot <===================================
plot1 <- create_scatter_plot(db_c, "37-38", "OR_37", "OR_38","SiU6atac", "Scr")
plot2 <- create_scatter_plot(db_c, "38-39", "OR_38", "OR_39","Scr", "ntc")
plot3 <- create_scatter_plot(db_c, "37-39", "OR_37", "OR_39","SiU6atac", "ntc")

plot4 <- create_scatter_plot(db_l, "46-47", "OR_46", "OR_47", "SiU6atac", "Scr")
plot5 <- create_scatter_plot(db_l, "46-48", "OR_46", "OR_48", "Scr", "ntc")
plot6 <- create_scatter_plot(db_l, "46-47", "OR_46", "OR_47", "SiU6atac", "Scr")

# Combine the plots together
combined_plots <- grid.arrange(plot1, plot2, plot3,                               
                               plot4, plot5, plot6, nrow = 2)
print(combined_plots)
# ========================== for the isoforms ==================================
create_scatter_plot_iso <- function(db, prefix, ortho_name1, ortho_name2,
                                    condition1, condition2){
  
  # Construct column names dynamically based on the prefix
  col_tpm_1 <- paste(prefix, "TPM_1", sep = " ")
  col_tpm_2 <- paste(prefix, "TPM_2", sep = " ")
  col_bin <- paste(prefix, "Bin", sep = " ")
  
  # Filter and select data for plotting
  filtered_db <- db %>%
    filter((
      !!rlang::sym(col_tpm_1) >= 0.5) | (!!rlang::sym(col_tpm_2) >= 0.5)) %>%
    select(all_of(c(col_tpm_1, col_tpm_2, col_bin)))
  
  # Create a scatter plot
  plot <- ggplot(filtered_db,
                 aes(x = log10(get(col_tpm_2)), y = log10(get(col_tpm_1)),
                     color = get(col_bin))) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = setNames(c("#c2d5f4", "#ff8c8c"),
                                         c(ortho_name1, ortho_name2))) +
    labs(title = paste("Scatter Plot of Log10 Transformed TPM Values (",
                       condition1,"vs", condition2, "), isoform"),
         x = paste("Log10(TPM_2) (", condition2, ")"),
         y = paste("Log10(TPM_1) (", condition1, ")"),
         color = "Isoform Bin") +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "white")) +
    theme(legend.position = "right")
  # Save the plot
  file_name <- paste("Log10_TPM_Scatter_Plot_", condition1,"vs",
                     condition2 , ".png", sep="")
  ggsave(file_name, plot, width = 8, height = 6)
  # Print the plot
  print(plot)
}
# ===========================> Usage example <==================================
plot1i <- create_scatter_plot_iso(db_c_iso, "37-38", "OR_37", "OR_38",
                                  "SiU6atac", "Scr")
plot2i <- create_scatter_plot_iso(db_c_iso, "38-39", "OR_38", "OR_39",
                                  "Scr", "ntc")
plot3i <- create_scatter_plot_iso(db_c_iso, "37-39", "OR_37", "OR_39",
                                  "SiU6atac", "ntc")
plot4i <- create_scatter_plot_iso(db_l_iso, "46-47", "OR_46", "OR_47",
                                  "SiU6atac", "Scr")
plot5i <- create_scatter_plot_iso(db_l_iso, "46-48", "OR_46", "OR_48",
                                  "Scr", "ntc")
plot6i <- create_scatter_plot_iso(db_l_iso, "46-47", "OR_46", "OR_47",
                                  "SiU6atac", "Scr")
combined_plots_iso <- grid.arrange(plot1i, plot2i, plot3i,
                                   plot4i, plot5i, plot6i, nrow = 2)
print(combined_plots_iso)
