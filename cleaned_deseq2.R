library(tidyverse)
library(miscTools)
library(readxl)
library(clusterProfiler)
library(DESeq2)
library('org.Mm.eg.db')
library("org.Hs.eg.db")
library('enrichplot')
library('ggplot2')
library('pheatmap')
library('plotCounts')
library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
#remotes::install_github("RubD/Giotto@cless")
library(Giotto)
library(dplyr)
library(msigdbr)
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

# all together 
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


# ==============================================================================
# Step 1: Select relevant columns
colData <- annot_data %>%
  dplyr::select(SlideName, SegmentDisplayName, U6ATAC_RNAscopeClass, SegmentLabel) %>%
  distinct()  # Ensure unique rows if needed

# Step 2: Create a combined condition factor (if needed)
colData$condition <- with(colData, paste(U6ATAC_RNAscopeClass, SegmentLabel, sep="_"))

# Step 3: Set rownames to match count data columns (assuming SegmentDisplayName matches)
rownames(colData) <- colData$SegmentDisplayName

# Print the head of colData to verify
print(head(colData))

countData <- target_count
countData <- as.data.frame(countData)
rownames(countData) <- countData[, 1]
countData <- countData[,-1 ]
correct_order <- colData$SegmentDisplayName
countData <- countData[, correct_order]
all(colnames(countData) == correct_order)  # This should return TRUE

countData <- as.matrix(countData)
countData <- round(countData)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# Running the DESeq analysis
dds <- DESeq(dds)
dds_vst<- vst(dds, blind=TRUE)
pca_plot <- plotPCA(dds_vst, intgroup = "condition", returnData = FALSE)
print(pca_plot + ggtitle("PCA of VST Transformed Data"))
write.csv(as.data.frame(results_dds), file = "u6atac_deseq2_results.csv")


dds_tumour <- results(dds, contrast=c("condition", "high_tumour","low_tumour"))
dds_tme <- results(dds, contrast=c("condition", "high_TME", "low_TME"))
dds_lowatac <- results(dds, contrast=c("condition", "low_tumour","low_TME"))
dds_highatac <- results(dds, contrast=c("condition", "high_tumour", "high_TME"))
# ------------------------------------------------------------------------------

results_dds <- dds_lowatac

results_df <- as.data.frame(results_dds)
results_df$log2FoldChange <- results_df$log2FoldChange
results_df$negLog10Pvalue <- -log10(results_df$pvalue)
# Add a column to mark E2F family genes
e2f_genes <- c("E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8")  
#check
rownames(results_df)
print(rownames(results_df)[rownames(results_df) %in% e2f_genes])
# Adding the highlight column correctly
results_df$highlight <- ifelse(rownames(results_df) %in% e2f_genes, "E2F Family", "Other")
# Filter to get only E2F Family genes
e2f_results <- results_df[results_df$highlight == "E2F Family", ]
# Print the results for E2F Family genes
print(e2f_results)
results_df$color_group <- ifelse(results_df$highlight == "E2F Family" & results_df$pvalue < 0.05, 
                                 "E2F Family Significant", 
                                 ifelse(results_df$highlight == "E2F Family", "E2F Family Not Significant", "Other"))


# Check the first few rows to confirm 'highlight' column exists
head(results_df)
results_df["E2F7", ]
# Create the basic volcano plot with adjusted text labels
volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = negLog10Pvalue)) +
  geom_point(data = subset(results_df, color_group == "Other"),
             aes(color = color_group), alpha = 0.5, size = 3) +
  geom_point(data = subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")),
             aes(color = color_group), size = 3.5) +
  scale_color_manual(values = c("Other" = "grey", 
                                "E2F Family Significant" = "red", 
                                "E2F Family Not Significant" = "#b48080")) +
  labs(title = "Low U6atac, tumour vs TME",
       x = "Log2 Fold Change", y = "-Log10 p-value",
       color = "Gene Group") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text_repel(data = subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")),
                  aes(label = rownames(subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")))),
                  size = 3, color = "black", box.padding = 0.35, point.padding = 0.5,
                  min.segment.length = 0, max.overlaps = Inf)
ggsave("Low U6atac, tumour vs TME volcano_plot_v2.pdf", volcano_plot, width = 10, height = 8, units = "in")

# Print the plot
print(volcano_plot)
dev.off()

# ==============================================================================
# TME
results_dds <- dds_tme

# Extract the log2 fold changes and -log10 of the p-values for plotting
results_df <- as.data.frame(results_dds)
results_df$log2FoldChange <- results_df$log2FoldChange  # correct any issues with log transformations
results_df$negLog10Pvalue <- -log10(results_df$pvalue)

# Add a column to mark E2F family genes
e2f_genes <- c("E2F1", "E2F2", "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8")  # Adjust list as necessary
rownames(results_df)
# Check if E2F genes are in your dataset
print(rownames(results_df)[rownames(results_df) %in% e2f_genes])
# Adding the highlight column correctly
results_df$highlight <- ifelse(rownames(results_df) %in% e2f_genes, "E2F Family", "Other")
# Add a new column to distinguish significant E2F genes from non-significant ones
results_df$color_group <- ifelse(results_df$highlight == "E2F Family" & results_df$pvalue < 0.05, 
                                 "E2F Family Significant", 
                                 ifelse(results_df$highlight == "E2F Family", "E2F Family Not Significant", "Other"))

# View the head to ensure the column is added correctly
head(results_df)

# Filter to get only E2F Family genes
e2f_results <- results_df[results_df$highlight == "E2F Family", ]
# Print the results for E2F Family genes
print(e2f_results)

# Check the first few rows to confirm 'highlight' column exists
head(results_df)
results_df["E2F7", ]


# Continue using the `color_group` as defined previously:
# "E2F Family Significant", "E2F Family Not Significant", "Other"

# Create the volcano plot
volcano_plot <- ggplot(results_df, aes(x = log2FoldChange, y = negLog10Pvalue)) +
  geom_point(data = subset(results_df, color_group == "Other"),
             aes(color = color_group), alpha = 0.5, size = 3) +
  geom_point(data = subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")),
             aes(color = color_group), size = 3.5) +
  scale_color_manual(values = c("Other" = "grey", 
                                "E2F Family Significant" = "red", 
                                "E2F Family Not Significant" = "#b48080")) +
  labs(title = "TME environment, high vs low U6atac",
       x = "Log2 Fold Change", y = "-Log10 p-value",
       color = "Gene Group") +
  theme_minimal() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_text_repel(data = subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")),
                  aes(label = rownames(subset(results_df, color_group %in% c("E2F Family Significant", "E2F Family Not Significant")))),
                  size = 3, color = "black", box.padding = 0.35, point.padding = 0.5,
                  min.segment.length = 0, max.overlaps = Inf)
ggsave("TME highU6atac volcano_plot v2.pdf", volcano_plot, width = 10, height = 8, units = "in")
# ========================= Top 10 down and up regu ===========================
# Order results by p-value and by absolute log2FoldChange to find the most significant changes
results_df <- as.data.frame(dds_tumour)
results_df$gene <- rownames(results_df)  # Ensure gene names are available as a column for easy referencing
# Assuming results_df is already created and contains the 'negLog10Pvalue' column

# Add gene identifiers from rownames as a column if not already present
results_df$gene <- rownames(results_df)

# Select the top 10 genes with the highest -log10(pvalue)
top_significant_genes <- results_df %>%
  arrange(desc(pvalue)) %>%
  head(10)
top_significant_genes


# ================================== GO! =======================================
# Assuming 'results_df' contains your DESeq2 results\
results_df <- dds_tumour
significant_genes <- rownames(subset(results_df, padj < 0.05))
# Convert gene symbols to Entrez IDs
entrez_ids <- bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
#n In bitr(significant_genes, fromType = "SYMBOL", toType = "ENTREZID",  :
# 0.65% of input gene IDs are fail to map...
# Perform GO enrichment analysis
ego <- enrichGO(gene          = entrez_ids$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",  # Biological Process; can also be "CC" for Cellular Component or "MF" for Molecular Function
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)

# Viewing the results
print(summary(ego))
barplot(ego, showCategory = 20) + 
  theme(
    text = element_text(size = 10),  # Adjusts global text size
    axis.text.x = element_text(size = 8),  # Adjusts x-axis text size
    axis.text.y = element_text(size = 8),  # Adjusts y-axis text size
    legend.text = element_text(size = 8)  # Adjusts legend text size
  )
# Dot plot for GO terms
dotplot(ego, showCategory = 20) +
  theme(
    text = element_text(size = 10),  # Adjusts global text size
    axis.text.x = element_text(size = 8),  # Adjusts x-axis text size
    axis.text.y = element_text(size = 8),  # Adjusts y-axis text size
    legend.text = element_text(size = 8)  # Adjusts legend text size
  )
# You can adjust the number of categories
# Create an enrichment map
cnetplot(ego, foldChange = entrez_ids$foldChange)  # Assuming fold changes are available; omit if not
# Heatmap of GO terms (if multiple GO analyses available)
# Not directly applicable here without multiple ego objects but included for completeness
heatmap(ego)


# ---------------- proliferation----------------------
library(dplyr)

# Assuming 'ego' is your enrichGO object
ego_results <- as.data.frame(ego)
ego_results$log2FoldChange <- ego_results$log2FoldChange  # correct any issues with log transformations
ego_results$negLog10Pvalue <- -log10(ego_results$pvalue)
proliferation_terms <- ego_results %>%
  filter(grepl("proliferation", Description, ignore.case = TRUE))  # Search for 'proliferation' in term descriptions

# View the proliferation related GO terms
print(proliferation_terms)
# Generate the dot plot and highlight the specific term
dotplot(ego, showCategory = 30) +  # Adjust showCategory as needed to display more terms
  geom_hline(data = subset(proliferation_terms, Description == "cellular proliferation"), 
             aes(yintercept = negLog10Pvalue), linetype = "dashed", color = "red", size = 1.5) +
  theme(
    axis.text.x = element_text(size = 6),
    axis.text.y = element_text(size = 6, color = ifelse(ego_results$ID %in% proliferation_terms$ID, "red", "black")),
    legend.text = element_text(size = 6)
  )


# ====================== Get Hallmark Genes ====================================
# Fetch Hallmark gene sets
#Handling failed to map genes:
library(clusterProfiler)

# Convert gene symbols to Entrez IDs
entrez_conversion <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# If multiple Entrez IDs are returned for a single gene symbol, keep only the first one
entrez_conversion <- entrez_conversion[!duplicated(entrez_conversion$SYMBOL), ]

# Merge the conversion table with the results dataframe
# This will add the Entrez IDs as a new column to your results_df
results_df <- merge(results_df, entrez_conversion, by.x = "row.names", by.y = "SYMBOL", all.x = TRUE)

# Replace NA in Entrez IDs with a placeholder or remove these rows
results_df$ENTREZID[is.na(results_df$ENTREZID)] <- "NA"
# Check how the merge went
print(head(results_df))

msig_data <- msigdbr(species = "Homo sapiens", category = "H")
# Convert gene symbols to Entrez IDs
gene_symbols <- rownames(results_df)
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Merge and clean data
results_df$entrez <- entrez_ids$ENTREZID
results_df <- na.omit(results_df)  # Remove rows with NA Entrez IDs

# Prepare gene list: Entrez IDs as names and statistics as values, sorted by significance
gene_list <- with(results_df, setNames(log2FoldChange, entrez))
gene_list <- sort(gene_list, decreasing = TRUE)  # Assuming higher values are more significant


