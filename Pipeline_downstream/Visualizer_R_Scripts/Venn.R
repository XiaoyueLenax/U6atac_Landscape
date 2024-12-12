library(VennDiagram)
library(RColorBrewer)
#  =============================================================================
#  Author: Xiaoyue Deng                       Contact: xiaoyued.deng@gmail.com
#  =============================================================================

# Get data: ====================================================================
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
#===============================================================================
#For LnCap
data4 <- read_and_rename_gene("IsoDE_raw_LnCap/46-47-DECalls-FC_gene.txt")
data5 <- read_and_rename_gene("IsoDE_raw_LnCap/46-48-DECalls-FC_gene.txt")
data6 <- read_and_rename_gene("IsoDE_raw_LnCap/47-48-DECalls-FC_gene.txt")
# now do the same for isoforms
data4a <- read_and_rename_isoform(
  "IsoDE_raw_LnCap/46-47-DECalls-FC_isoform.txt")
data5a <- read_and_rename_isoform(
  "IsoDE_raw_LnCap/46-48-DECalls-FC_isoform.txt")
data6a <- read_and_rename_isoform(
  "IsoDE_raw_LnCap/47-48-DECalls-FC_isoform.txt")

# Define a function that filter based on expression threshold-------------------
filter_expthr <- function(db, prefix){
  # Construct column names dynamically based on the prefix
  col_tpm_1 <- paste(prefix, "TPM_1", sep = " ")
  col_tpm_2 <- paste(prefix, "TPM_2", sep = " ")
filtered_db <- db %>%
  filter((.[[col_tpm_1]] >= 0.5 | .[[col_tpm_2]] >= 0.5))
}
# ---------on data----------------------
# These are the ones with elevated expression in both. 
f1 <- filter_expthr(data1, "37-38")
f2 <- filter_expthr(data2, "37-39")
f3 <- filter_expthr(data3, "38-39")
f4 <- filter_expthr(data4, "46-47")
f5 <- filter_expthr(data5, "46-48")
f6 <- filter_expthr(data6, "47-48")
# To find in one or another?
# Do the diagram ===============================================================
col <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(f1$Gene, f2$Gene, f3$Gene),
  main = "Venn Diagram of IsoDE binning outputs (Testing)",
  main.col = "#9600ff",
  category.names = c("SiU6atac - Scr" , "SiU6atac - ntc " , "Scr - ntc"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  fill = col
)

# Try a big one:
col <- brewer.pal(6, "Pastel2")
venn.diagram(
  x = list(f1$Gene, f2$Gene, f3$Gene, f4$Gene,f5$Gene, f6$Gene),
  main = "Venn Diagram of IsoDE binning outputs (Testing)",
  main.col = "#9600ff",
  category.names = c( "C42 SiU6atac - Scr" , "C42 SiU6atac - ntc " , "C42Scr - ntc",
                      "LnCap U6 - Scr", "LnCap U6 - ntc", "LnCap Scr - ntc"),
  filename = 'all_venn_diagramm.png',
  output=TRUE,
  fill = col
)

col <- brewer.pal(3, "Pastel2")
venn.diagram(
  x = list(f1$Gene, f4$Gene),
  main = "Venn Diagram of IsoDE binning outputs (Testing)",
  main.col = "#9600ff",
  category.names = c( "C42" , "LnCap"),
  filename = 'all_venn_diagramm.png',
  output=TRUE,
  fill = c("#eeaf61","#6a0d83")
)
