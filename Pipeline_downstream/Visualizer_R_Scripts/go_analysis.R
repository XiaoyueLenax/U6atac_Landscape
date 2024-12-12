library(readxl)
library(openxlsx)
library(dplyr)
library(readr)
library(ggplot2)
library(gridExtra)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(KEGGREST)
library(ggrepel)
library(pathfindR)

# load the dataset
db_c <- read_excel("Binning_stats_C42_Genes.xlsx")
db_l <- read_excel("Binning_stats_LnCap_Genes.xlsx")
db_l_iso <- read_excel("Binning_stats_LnCap_isoforms.xlsx")
db_c_iso <- read_excel("Binning_stats_C42_isoforms.xlsx")
project_title <- 'C42'

#==========================
#------ KEGG Plot ---------
#==========================
get_filtered_gene <- function(db, prefix, threshold_a, threshold_b) {
  # Construct column names dynamically based on the prefix
  col_tpm_1 <- paste(prefix, "TPM_1", sep = " ")
  col_tpm_2 <- paste(prefix, "TPM_2", sep = " ")
  col_bin <- paste(prefix, "Bin", sep = " ")
  
  # Filter and select data for plotting, including the gene ID column
  filtered_db <- db %>%
    dplyr::filter((!!rlang::sym(col_tpm_1) >= threshold_a) | (!!rlang::sym(col_tpm_2) >= threshold_b)) %>%
    dplyr::select(1, all_of(c(col_tpm_1, col_tpm_2, col_bin)))  # Include the first column (gene ID) in the selection
  
  return(filtered_db)
}

# Example usage:
sample <- "C42_37_38"
date <-Sys.Date()
gene_list <- get_filtered_gene(db_c, '37-38', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_c, '37-39', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_c, '38-39', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l, '46-47', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l, '46-48', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l, '47-48', 0.5, 0.5)

#gene_list <- get_filtered_gene(db_c_iso, '37-38', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_c_iso, '37-39', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_c_iso, '38-39', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l_iso, '46-47', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l_iso, '46-48', 0.5, 0.5)
#gene_list <- get_filtered_gene(db_l_iso, '47-48', 0.5, 0.5)



#convery with entriz ID
entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list$Gene...1, column = "ENTREZID", keytype = "SYMBOL")
collapsed_entrez_ids <- sapply(entrez_ids, function(x) ifelse(length(x) > 0, x[1], NA))
collapsed_entrez_ids

# Perform GO enrichment analysis
kegg_enrichment <- enrichKEGG(gene = collapsed_entrez_ids,
                              organism = "hsa",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.05)
print(kegg_enrichment)
#enrichment_chart(kegg_enrichment)
#Lower q-values indicate greater statistical significance and are often associated with enriched pathways or terms. 
top_enriched <- kegg_enrichment[kegg_enrichment$qvalue < 0.05, ]
top_enriched <- top_enriched[order(top_enriched$qvalue), ]
head(top_enriched)
tail(top_enriched)
# browse pathway
browseKEGG(kegg_enrichment,"hsa04141")
browseKEGG(kegg_enrichment,"hsa04140")
browseKEGG(kegg_enrichment,"hsa04110")
browseKEGG(kegg_enrichment,"hsa04144")



png(paste0("top KEGG pathways ", date, "_", sample, ".png"), width = 1800, height = 3000)
#colors <- scales::gradient_n_pal(c("#f7cac9", "#92a8d1"))(top_enriched$qvalue)
# Plot the top enriched pathways
ggplot(top_enriched, aes(x = -log10(qvalue), y = Description)) +
  geom_point(size = 6, color = "#abb1cf") +
  geom_segment(aes(xend = 0, yend = Description), color = "#c5b9cd", lwd=3) +
  labs(x = "-log10(Q-value)", y = "Pathway", title = "Top Enriched KEGG Pathways") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 20)) +  # Adjust y-axis label size
  xlim(0, max(-log10(top_enriched$qvalue) + 1))  # Set x-axis limit
dev.off()

#-------- another visualizer? ------------
library("pathview")
hsa04110 <- pathview(gene.data  = collapsed_entrez_ids,
                     pathway.id = "hsa04110",
                     species    = "hsa")
hsa04144 <- pathview(gene.data  = collapsed_entrez_ids,
                     pathway.id = "hsa04144",
                     species    = "hsa")
hsa05235 <- pathview(gene.data  = collapsed_entrez_ids,
                     pathway.id = "hsa05235",
                     species    = "hsa")
# --------------------- END OF KEGG -------------------
# ====================================================
# ---- Enrich GO plots -------------------------------
# ===================================================
# since the gene one does not have emsembl ID, use SYMBOL 

go_enrich_result <- enrichGO(
  gene = collapsed_entrez_ids,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process ontology
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
go_enrich_result
enriched_terms <- go_enrich_result@result
quota = 50
significant_terms <- enriched_terms[enriched_terms$p.adjust < 0.05, ]
png(paste0("Top", quota, "Enriched Terms  ", project_title, date,".png"), width = 900, height = 1500)
mutate(go_enrich_result, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",showCategory=quota)
dev.off()


# dot plot 
png(paste0("Top", quota, "Enriched Terms Dot Plot", project_title, date,".png"), width = 1500, height = 1000)
dotplot(go_enrich_result, showCategory=30) + ggtitle("dotplot for",project_title)
dev.off()


#attempt at network
#p1 <- cnetplot(go_enrich_result, foldChange=db_c)
