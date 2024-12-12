library(tidyverse)
library(miscTools)
library(readxl)
library(clusterProfiler)
library(DESeq2)
library("org.Hs.eg.db")
library('enrichplot')
library('ggplot2')
library('pheatmap')
library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
library(Giotto)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(openxlsx)
library(msigdbr)
#Getting All the libraries prepared
library(tidyverse)
library(ggplot2)
library(DOSE)
library(ggupset)
library(DESeq2)
library(ggrepel)
library(EnhancedVolcano)
library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
library(msigdbr)
library(ReactomePA)
library(EnhancedVolcano)
library(umap)

# sanity check for signature markers to ensure data structure/ 
# gsva analysis 
# UMAP
# transcription factor analysis 
# GSEA analysis from cancer genomics
#===============================================================================
# Author: Xiaoyue Deng                  Contact: xiaoyue.deng@students.unibe.ch
#============================== General Data prep ==============================
annot_data <- read_excel("1_4-12-2024/prostate_rawdata_100424_annon.xlsx",
                         sheet = 1)
target_properties <- read_excel("1_4-12-2024/prostate_rawdata_100424_annon.xlsx",
                                sheet = 2)
target_count <- read_excel("1_4-12-2024/prostate_rawdata_100424_annon.xlsx",
                           sheet = 3)
bio_probe <- read_excel("1_4-12-2024/prostate_rawdata_100424_annon.xlsx", sheet = 4)
# warning Warning message:                                                                                                                                                         
#Expecting logical in I11505 / R11505C9: got 'Failed Grubbs test and excluded from all AOIs, AOI proportion is: 3.061224E-01'
db_sum <- read_excel("1_4-12-2024/prostate_rawdata_100424_annon.xlsx", sheet = 5)
annot_data <- annot_data %>% mutate("pAlignedReads"=AlignedReads/RawReads*100)
annot_data <- annot_data %>% mutate("ExpressionFilteringThreshold..Human.NGS.Whole.Transcriptome.Atlas.RNA_1.0."=
                                      ifelse("LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0"<=2, 2,
                                             ifelse("LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0">2,
                                                    "LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0", NA)))

# -------------------------------- Further Trimming ----------------------------
annot_data_filtered <- annot_data[annot_data$`ExpressionFilteringThreshold (Human NGS Whole Transcriptome Atlas RNA_1.0)` > 2,]
annot_data <- annot_data_filtered
# Basic: DESeq =================================================================
colData <- annot_data %>%
  dplyr::select(SlideName, SegmentDisplayName, U6ATAC_RNAscopeClass, SegmentLabel)
colData$condition <- with(colData, paste(U6ATAC_RNAscopeClass, SegmentLabel, sep="_"))
rownames(colData) <- colData$SegmentDisplayName

print(head(colData))

countData <- target_count
countData <- as.data.frame(countData)
rownames(countData) <- countData[, 1]
countData <- countData[,-1 ]
correct_order <- colData$SegmentDisplayName
countData <- countData[, correct_order]
all(colnames(countData) == correct_order)  # This should return TRUE

# makr matrix and round them. 
countData <- as.matrix(countData)
countData <- round(countData)
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
# Running the DESeq2 analysis
dds <- DESeq(dds)
dds_vst<- vst(dds, blind=TRUE)
pca_plot <- plotPCA(dds_vst, intgroup = "condition", returnData = FALSE)
print(pca_plot + ggtitle("PCA of VST Transformed Data, strict threshold"))

# slice our conditions of interest _!!!!!!! change contrast group, only U6atac.

# filtering step - cook's cutoff filtering, a parameter to set in the results section.  
dds_tumour <- results(dds, contrast=c("condition", "high_tumour","low_tumour"))
# high pri vs low pri
dds_tme <- results(dds, contrast=c("condition", "high_TME", "low_TME"))
#dds_lowatac <- results(dds, contrast=c("condition", "low_tumour","low_TME"))
#dds_highatac <- results(dds, contrast=c("condition", "high_tumour", "high_TME"))
#=====================

# Slice out a gene

#====================
# Shows a specific gene count
gene_name <- "PTPRC"
counts <- plotCounts(dds, gene = gene_name, intgroup = "condition", returnData = TRUE)
counts <- tibble::rownames_to_column(counts, "SegmentDisplayName")
pos <- position_jitter(width = 0.3, seed = 2)

# Plot
counts %>%
  ggplot(aes(x = condition, y = count, fill = condition)) +
  geom_boxplot(outlier.color = NA) +
  geom_jitter(position = pos) +
  theme_bw() +
  theme(axis.title.x.bottom = element_blank(), legend.position = "none") +
  labs(title = paste("Counts of", gene_name, "by condition, very soft threshold"))  # Adjusted title


# ==============================================================================
#           CANCER GENOMICS 
# ==============================================================================
results_df<- as.data.frame(dds_tumour)
tag <- "tumour"
#tag <- "TME"
# volcano - wip
EnhancedVolcano(results_df,
                lab = rownames(results_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = paste(tag, " high vs low U6atac"),
                legendPosition = 'bottom')

# GSEA analysis ==== lengthy conversion ========================================
results_df_gsea <- results_df %>% rownames_to_column(var="SYMBOL")
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=results_df_gsea$SYMBOL, 
                                    columns=("ENTREZID"), 
                                    keytype="SYMBOL") 
ens2symbol <- dplyr::as_tibble(ens2symbol)
res <- dplyr::inner_join(ens2symbol, results_df_gsea, by="SYMBOL")
res2 <- res %>% 
  dplyr::select(ENTREZID, log2FoldChange) %>% 
  na.omit() %>% 
  dplyr::distinct() %>% 
  group_by(ENTREZID) %>%  ## this averages the log2fc across duplicate genes
  dplyr::summarize(stat=mean(log2FoldChange)) # stat variable is now equivalent to log2fc

res2 <- as.data.frame(res2)
geneList = res2[,2] # numerical variable, gene counts, log2fc
names(geneList) = as.character(res2[,1]) # character variable for gene names
geneList = sort(geneList, decreasing = TRUE) # organizes genelist descending log2fc

## THIS IS OUR GENESETS AND PATHWAYS TO ASSESS
H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# perform GSEA
em2 <- GSEA(geneList, TERM2GENE = H)

## visualise GSEA DOTPLOT
category_num = 30
dot1 <- dotplot(em2, showCategory=category_num, split=".sign", font.size=8, size = NULL, label_format = 90) + 
  ggtitle(paste("top",category_num, " hallmarks of", tag, "high vs low U6atac")) +
  facet_grid(~.sign)
dot1 + 
  theme(plot.title.position = "plot")


filename <- paste("top", category_num, "hallmarks of", tag, "high vs low U6atac", ".png", sep = "_")
ggsave(filename, plot = dot1, device = "png")

## visualise specific pathways 
em2@result[["ID"]]
gseaplot2(em2, geneSetID = 8, title = em2$Description[8]) #Works
temp_gsea <- em2@result
## OVER REPRESENTATION ANALYSIS AND VISUALISATION
data(geneList, package="DOSE")
gene <- names(geneList)
em <- enricher(gene, TERM2GENE=H)
## visualise ORA
upsetplot(em) #Works
# ==============================================================================
#           TME
# ==============================================================================
results_df<- as.data.frame(dds_tme)
tag <- "TME"
#tag <- "TME"
# volcano - wip
EnhancedVolcano(results_df,
                lab = rownames(results_df),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = paste(tag, " high vs low U6atac"),
                legendPosition = 'bottom')

# GSEA analysis ==== lengthy conversion ========================================
results_df_gsea <- results_df %>% rownames_to_column(var="SYMBOL")
ens2symbol <- AnnotationDbi::select(org.Hs.eg.db,
                                    key=results_df_gsea$SYMBOL, 
                                    columns=("ENTREZID"), 
                                    keytype="SYMBOL") 
ens2symbol <- dplyr::as_tibble(ens2symbol)
res <- dplyr::inner_join(ens2symbol, results_df_gsea, by="SYMBOL")
res2 <- res %>% 
  dplyr::select(ENTREZID, log2FoldChange) %>% 
  na.omit() %>% 
  dplyr::distinct() %>% 
  group_by(ENTREZID) %>%  ## this averages the log2fc across duplicate genes
  dplyr::summarize(stat=mean(log2FoldChange)) # stat variable is now equivalent to log2fc

res2 <- as.data.frame(res2)
geneList = res2[,2] # numerical variable, gene counts, log2fc
names(geneList) = as.character(res2[,1]) # character variable for gene names
geneList = sort(geneList, decreasing = TRUE) # organizes genelist descending log2fc

## THIS IS OUR GENESETS AND PATHWAYS TO ASSESS
H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

# perform GSEA
em2 <- GSEA(geneList, TERM2GENE = H)

## visualise GSEA DOTPLOT
category_num = 30
dot1 <- dotplot(em2, showCategory=category_num, split=".sign", font.size=8, size = NULL, label_format = 90) + 
  ggtitle(paste("top",category_num, " hallmarks of", tag, "high vs low U6atac")) +
  facet_grid(~.sign)
dot1 + 
  theme(plot.title.position = "plot")


filename <- paste("top", category_num, "hallmarks of", tag, "high vs low U6atac", ".png", sep = "_")
ggsave(filename, plot = dot1, device = "png")

## visualise specific pathways 
em2@result[["ID"]]
gseaplot2(em2, geneSetID = 8, title = em2$Description[8]) #Works
temp_gsea <- em2@result
## OVER REPRESENTATION ANALYSIS AND VISUALISATION
data(geneList, package="DOSE")
gene <- names(geneList)
em <- enricher(gene, TERM2GENE=H)
## visualise ORA
upsetplot(em) #Works
# ===================== between group comparisons? =============================
#                          AR Score Analysis
# ==============================================================================
ar_genes <- read_excel("AR_NEPC_scores.xlsx",sheet=2)
ar_genes <- ar_genes$Trimmed_AR
ar_genes
# Remove NAs from ar_genes
ar_genes <- na.omit(ar_genes)
nrow(umap_data)

sum(rownames(results_df) %in% ar_genes)
ar_gene_expr <- results_df[rownames(results_df) %in% ar_genes, ]
ar_z_scores <- apply(ar_gene_expr, 1, function(x) (x - mean(x)) / sd(x))
ar_z_scores
ar_gene_expr$z_score <- ar_z_scores
# calculate with rawreads for normalization to get AR score !!!!
# figures spider plot, normalied gsva (averaged across high and low U6atac) - rank within each group the enrichment score,
# assign values to pathway 
#spatial-omics overlay
#nanostring superimpose analysis 

# spatial deconvolution. 
# library of datasets that defines this cell type. 
# expect to see enrich keratin markers. CD45


# set up meeting time with Hannah. 
# Monday 10:30.