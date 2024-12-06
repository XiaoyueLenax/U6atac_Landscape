library(pheatmap)
library(dplyr)
library(SpatialDecon)
browseVignettes("SpatialDecon")
library(SeuratObject)
library(NanoStringGeoMxTools)
library(DESeq2)
#===============================================================================
#                           CONTINUING FROM DESEQ2
#===============================================================================
rld <- rlog(dds, blind = TRUE)
vst <- vst(dds, blind = TRUE)
rlog_counts <- assay(rld)
vst_counts <- assay(vst)
class(vst_counts)
dim(rlog_counts)
nrow(rlog_counts)
head(rlog_counts)
#checking for negprobe
#=============================== investigation into negprob
bg_pro = derive_GeoMx_background(norm=rlog_counts, probepool = rep(1, nrow(rlog_counts)), negnames = "NegProbe-WTX")
#bg_pro = derive_GeoMx_background(norm=vst, probepool = rep(1, nrow(vst_counts)), negnames = "NegProbe-WTX")
#vst
# Check the structure of prostate_data
#-------- loading prostate data
load("Prostate_Henry.RData")
write.csv(x = profile_matrix, file = "Prostate_matrix.csv", 
          row.names = TRUE, quote = FALSE)
profile_matrix <- as.matrix(profile_matrix)
str(profile_matrix)

# Assuming annot_data is already loaded and contains the necessary information

# Perform spatial deconvolution
decon_results <- spatialdecon(norm = rlog_counts,
                              bg = bg_pro,
                              X = profile_matrix,
                              align_genes = TRUE)

# Handle missing values (replace NA with row/column mean)
rlog_counts[is.na(rlog_counts)] <- rowMeans(rlog_counts, na.rm = TRUE)
rlog_counts[is.na(rlog_counts)] <- colMeans(rlog_counts, na.rm = TRUE)

# Extract U6ATAC_RNAscopeClass from annot_data
u6atac_class <- as.factor(annot_data$U6ATAC_RNAscopeClass)
head(u6atac_class)
head(annot_data$SegmentDisplayName)
# Create an annotation data frame
annotation_data <- data.frame(U6ATAC_RNAscopeClass = u6atac_class)
rownames(annotation_data) <- annot_data$SegmentDisplayName
# Calculate mean values for each row (cell type) in decon_results$beta
mean_values <- rowMeans(decon_results$beta, na.rm = TRUE)

# Order the rows by mean values in descending order
sorted_rows <- order(mean_values, decreasing = TRUE)
sorted_beta <- decon_results$beta[sorted_rows, ]  # Sort beta matrix by row means

# Sort the columns based on U6ATAC_RNAscopeClass
sorted_columns <- order(u6atac_class)
sorted_beta <- sorted_beta[, sorted_columns]  # Sort beta matrix by column classes

# Reorder the annotation data frame accordingly
sorted_annotation_data <- annotation_data[sorted_columns, , drop = FALSE]

# Define a color palette for the heatmap
color_palette <- colorRampPalette(c("#4f5bd5", "#d62976", "#feda75"))(50)

# Define the colors for the annotation
annotation_colors <- list(
  U6ATAC_RNAscopeClass = c("low" = "#3accbd", "high" = "#ff8b94")  # Define colors for annotations
)

# Plot the heatmap with sorted samples and rows
#pdf("heatmap_sorted_U6ATAC_Tumour_class.pdf", width = 20, height = 12)
pheatmap(sorted_beta,
           labels_row = rownames(decon_results$beta)[sorted_rows],  # Sorted row labels
           main = "Sorted Spatial Deconvolution Heatmap by U6ATAC Class",
           fontsize_row = 10,
           fontsize_col = 10,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           color = color_palette,
           border_color = "black",
           annotation_col = sorted_annotation_data,  # Add U6ATAC annotations
           annotation_colors = annotation_colors)  # Set colors for annotations
#dev.off()
   