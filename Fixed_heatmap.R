## Load necessary libraries for GSEA
library(clusterProfiler)
library(org.Hs.eg.db)     # Human gene annotation
library(msigdbr)          # MSigDB gene sets
library(enrichplot)       # For plotting enrichment results
library(dplyr)
library(ggplot2)
## Step 1: Prepare the Ranked Gene List from DESeq2 Results
# Assuming 'res' contains the DESeq2 results
res <- results(dds, contrast = c("U6ATAC_RNAscopeClass", "high", "low"))
res <- as.data.frame(res) %>% filter(!is.na(log2FoldChange))

# Create the gene list for GSEA
gene_list <- res$log2FoldChange
names(gene_list) <- rownames(res)  # Use gene symbols as names
gene_list <- sort(gene_list, decreasing = TRUE)

## Step 2: Load Hallmark Gene Sets and Prepare for GSEA
# Load MSigDB Hallmark gene sets for Homo sapiens
hallmark_gene_sets <- msigdbr(species = "Homo sapiens", category = "H")

# Convert `hallmark_gene_sets` into a `TERM2GENE` format data frame for GSEA
TERM2GENE <- hallmark_gene_sets %>% dplyr::select(gs_name, gene_symbol)

# Filter `gene_list` to include only genes found in `TERM2GENE`
gene_list <- gene_list[names(gene_list) %in% TERM2GENE$gene_symbol]

## Step 3: Run GSEA with Corrected Input
gsea_results <- GSEA(geneList = gene_list, TERM2GENE = TERM2GENE, pvalueCutoff = 0.05)

## Step 4: Plot Top 12 Enriched Pathways
# Extract and plot the top 12 enriched pathways
dotplot(gsea_results, showCategory = 10) +
  ggtitle("Top 10 Enriched Hallmark Pathways") +
  theme_minimal()
