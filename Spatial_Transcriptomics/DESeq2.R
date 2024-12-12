library(tidyverse)
library(miscTools)
library(readxl)
library(clusterProfiler)
library(DESeq2)
library('org.Mm.eg.db')
library('enrichplot')
library('ggplot2')
library('pheatmap')
library('plotCounts')
library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
#remotes::install_github("RubD/Giotto@cless")
library(Giotto)
library(dplyr)
#============================== General Data prep ==============================
annot_data <- read_excel("prostate_rawdata_100424_annon.xlsx",
                         sheet = 1)
target_properties <- read_excel("prostate_rawdata_100424_annon.xlsx",
                                sheet = 2)
target_count <- read_excel("prostate_rawdata_100424_annon.xlsx",
                           sheet = 3)
bio_probe <- read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 4)
# warning Warning message:
#Expecting logical in I11505 / R11505C9: got 'Failed Grubbs test and excluded from all AOIs, AOI proportion is: 3.061224E-01'
db_sum <- read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 5)
annot_data <- annot_data %>% mutate("pAlignedReads"=AlignedReads/RawReads*100)
#smart changing of the column names based on criteria?
annot_data <- annot_data %>% mutate("ExpressionFilteringThreshold..Human.NGS.Whole.Transcriptome.Atlas.RNA_1.0."=
                                      ifelse("LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0"<=2, 2,
                                             ifelse("LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0">2, "LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0", NA)))
#===============================================================================
# list of conditions
# tumor high u6atac vs tumor low u6atac
# TME high U6atac vs low U6atac
#  tumor vs tme high u6atac
#  tumor vs tme low u6atac
tumour_all <- annot_data[annot_data$SegmentLabel == "tumour", ]
tum_hi_u6 <- tumour_all[tumour_all$U6ATAC_RNAscopeClass == "high", ]
tum_low_u6 <- tumour_all[tumour_all$U6ATAC_RNAscopeClass == 'low', ]

sample_id <- tumour_all$SegmentDisplayName
siu6atac_conditions <- c(tumour_all$U6ATAC_RNAscopeClass)

reftable <- data.frame(siu6atac_conditions)
rownames(reftable) <- sample_id
high_u6_all <- annot_data[annot_data$U6ATAC_RNAscopeClass == "high",]
low_u6_all <- annot_data[annot_data$U6ATAC_RNAscopeClass == "low",]

hi_u6_names <- high_u6_all$SegmentDisplayName %>%
  intersect(colnames(target_count))
high_u6_counts_ <- target_count[, hi_u6_names]

hi_u6_condition <- data.frame(high_u6_all$SegmentLabel)
rownames(hi_u6_condition) <- high_u6_all$SegmentDisplayName
gene_names <- target_count[[1]]

hi_u6_counts <- cbind(TargetName = gene_names, target_count[, hi_u6_names])
head(gene_names)
head(high_u6_counts)

high_u6_counts <- round(high_u6_counts)
colnames(hi_u6_condition) <- "Conditions"
head(hi_u6_counts)

dds <- DESeqDataSetFromMatrix(countData = high_u6_counts,
                              colData = reftable,
                              design = ~ Conditions)

dds <- DESeq(dds)
dds$TargetName #Ok so, sill
res <- results(dds)
summary(res)
res_sig <- res[which(res$padj < 0.05), ]
#PCA plot
dds_vst<- vst(dds, blind=TRUE)
plotPCA(dds_vst, intgroup = 'Conditions', main = "High ")
write.csv(as.data.frame(res_sig), file = "high_u6atac_deseq2_results.csv")
#======================== Visualization (WIP) ==================================
dna_repair_genes <- c("BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "RAD51", "XRCC1",
                      "XRCC2", "XRCC3", "RAD51C", "RAD51D", "RAD52", "RAD54L",
                      "BRIP1", "BARD1", "NBN", "RAD50", "MRE11", "ATR", "ATRIP",
                      "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF",
                      "BRCC3", "FANCG", "FANCI", "UBE2T", "FANCL", "FANCM",
                      "PTEN", "MLH1", "MSH2", "MSH3", "MSH6", "PMS1",
                      "PMS2", "POLD1", "POLE", "EXO1", "RAD1", "RAD9A", "HUS1",
                      "RAD17", "TOPBP1", "MDC1")
hi_u6<-DESeq2::results(dds_output, contrast = c('conditions','high'))
# 1. Volcano plot
# Convert DESeqResults object to dataframe
res_df <- as.data.frame(res)
plot<-ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)",
       title = "Volcano Plot")
## other code from previous
library(ggplot2)
# Define color palette
point_colors <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0, "#E1B087", "gray")
# Create the plot
plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = point_colors, size = 3) +  # Using defined color palette and reducing point size
  geom_vline(xintercept = 0, col = "#4E4E4E", linetype = "dashed", size = 1) +  # Adjusting line settings
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)", title = "High U6atac") +  # Adjusting labels
  scale_color_manual(values = c("red", "gray")) +  # Manually set colors
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),  # Title formatting
    axis.title = element_text(size = 16, face = "bold"),  # Axis title formatting
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),  # Remove background grid
    axis.line.x = element_line(color = "black", size = 1),  # Solid x-axis line
    axis.line.y = element_blank(),  # Remove y-axis line
    axis.text = element_text(size = 16),  # Enlarge axis text
    axis.text.x = element_text(size = 16)
  )
# Define DNA repair genes
dna_repair_genes <- c("BRCA1", "BRCA2", "ATM", "CHEK2", "PALB2", "RAD51", "XRCC1",
                      "XRCC2", "XRCC3", "RAD51C", "RAD51D", "RAD52", "RAD54L",
                      "BRIP1", "BARD1", "NBN", "RAD50", "MRE11", "ATR", "ATRIP",
                      "FANCA", "FANCB", "FANCC", "FANCD2", "FANCE", "FANCF",
                      "BRCC3", "FANCG", "FANCI", "UBE2T", "FANCL", "FANCM",
                      "PTEN", "MLH1", "MSH2", "MSH3", "MSH6", "PMS1",
                      "PMS2", "POLD1", "POLE", "EXO1", "RAD1", "RAD9A", "HUS1",
                      "RAD17", "TOPBP1", "MDC1")
# Create a new column to indicate whether a gene is a DNA repair gene
res_df$is_dna_repair_gene <- ifelse(res_df$GeneSymbol %in% dna_repair_genes)
# Create the plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = is_dna_repair_gene)) +
  geom_point(alpha = 0.7, size = 3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  scale_color_manual(values = c("DNA Repair" = "blue", "Other" = "gray")) +  # Set color for DNA repair genes
  labs(x = "log2 Fold Change", y = "-log10(Adjusted p-value)", title = "Volcano Plot") +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    legend.position = "top"
  )
# Show the plot
print(volcano_plot)
plot  # Show the plot
# 2. MA plot
ggplot(res_df, aes(x = baseMean, y = log2FoldChange)) +
  geom_point() +
  labs(x = "Mean Expression", y = "log2 Fold Change",
       title = "MA Plot")
# 4. PCA plot
plotPCA(dds)
# 5. Clustered heatmap
# PCA plot
plotPCA(dds, intgroup = "Conditions")












