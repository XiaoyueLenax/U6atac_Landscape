## load libraries
library(tidyverse)
library(miscTools)
library(readxl)
library("org.Hs.eg.db")
library(devtools)  # if not installed: install.packages('devtools')
library(remotes)  # if not installed: install.packages('remotes')
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(openxlsx)
library(clusterProfiler)
library(DESeq2)
library('enrichplot')
library('ggplot2')
library("EnhancedVolcano")
library('pheatmap')
# This is the version that works for spatial decon, where we had to manually add back the
# Negprobe row.
# The original version of DSP does not work. 15th october 2024

setwd("~/Thesis/spatial_transcriptomics")
annot_data <- readxl::read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 1)
target_properties <- read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 2)
count_data <- read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 3)
bio_probe <- read_excel("prostate_rawdata_100424_annon.xlsx", sheet = 4)

tumour_annot <- annot_data[annot_data$SegmentLabel == "tumour", ]
tme_annot <- annot_data[annot_data$SegmentLabel == "TME", ]

#annot_data <- tumour_annot
annot_data <- tme_annot
annot_data <- annot_data %>% mutate("pAlignedReads" = AlignedReads / RawReads * 100)

## Ensure LOQ column exists
annot_data <- annot_data %>% mutate(
  "ExpressionFilteringThreshold..Human.NGS.Whole.Transcriptome.Atlas.RNA_1.0." =
    ifelse(
      "LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0" <= 2, 2,
      ifelse(
        "LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0" > 2,
        "LOQ (Human NGS Whole Transcriptome Atlas RNA_1.0", NA
      )
    )
)

## Filter out segments with <60% pAlignedReads
annot_data.1 <- annot_data %>% filter("pAlignedReads" >= 60)

## Select samples from count_data based on the filtered segments
count_data.1 <- count_data %>% dplyr::select(TargetName, annot_data.1$SegmentDisplayName)

################# Filtering steps for segments and targets #################

#### SEGMENT FILTERING
xl_metadata <- annot_data.1
xl_counts <- count_data.1

xl_annot <- xl_metadata %>%
  dplyr::select(
    SegmentDisplayName,
    "ExpressionFilteringThreshold (Human NGS Whole Transcriptome Atlas RNA_1.0)",
    SegmentLabel,
    pAlignedReads
  )

xl_counts <- xl_counts %>% column_to_rownames("TargetName")
xl_counts <- xl_counts[, match(xl_annot$"SegmentDisplayName", colnames(xl_counts))]

## Exclude samples with less than XX% of genes above LOQ
LOQ_sample_thresh <- 3
LOQ_sample_n <- round(nrow(xl_counts) / 100 * LOQ_sample_thresh)

counts_LOQ_bin <- t(t(xl_counts) > xl_annot$"ExpressionFilteringThreshold (Human NGS Whole Transcriptome Atlas RNA_1.0)")
counts_LOQ_sum_sample <- unname(apply(counts_LOQ_bin, 2, sum))

filtered_samples_1 <- xl_annot$SegmentDisplayName[counts_LOQ_sum_sample < LOQ_sample_n]

# Update the data frames
counts_LOQ_bin <- counts_LOQ_bin[, !(colnames(counts_LOQ_bin) %in% filtered_samples_1)]
xl_annot <- xl_annot[!(xl_annot$SegmentDisplayName %in% filtered_samples_1), ]
xl_counts <- xl_counts[, !(colnames(xl_counts) %in% filtered_samples_1)]

##### TARGET FILTERING
LOQ_gene_thresh <- 10
LOQ_gene_n <- round(ncol(xl_counts) / 100 * LOQ_gene_thresh)

counts_LOQ_sum_gene <- unname(apply(counts_LOQ_bin, 1, sum))
filtered_genes <- rownames(counts_LOQ_bin)[counts_LOQ_sum_gene < LOQ_gene_n]

# Update the data frames
counts_LOQ_bin <- counts_LOQ_bin[!(rownames(counts_LOQ_bin) %in% filtered_genes), ]
xl_counts <- xl_counts[!(rownames(xl_counts) %in% filtered_genes), ]

### ADD BACK NEGPROBE ROW ###
# Extract the NegProbe-WTX row from the original data (count_data.1)
negprobe_row <- count_data.1 %>%
  filter(TargetName == "NegProbe-WTX") %>%
  column_to_rownames("TargetName")

# Ensure column order matches between negprobe_row and xl_counts
negprobe_row <- negprobe_row[, colnames(xl_counts)]

# Add the NegProbe-WTX row back to xl_counts
xl_counts <- rbind(xl_counts, negprobe_row)

# Check if the row was added correctly
grep("NegProbe-WTX", rownames(xl_counts))  # This should return the row index

# Add TargetName back as a column for later processing
xl_counts <- tibble::rownames_to_column(xl_counts, "TargetName")

# Continue with the rest of the dataset preparation
xl_annot$U6ATAC_RNAscopeClass <- annot_data.1$U6ATAC_RNAscopeClass[match(xl_annot$SegmentDisplayName, annot_data.1$SegmentDisplayName)]
#xl_annot$PatientNum <- annot_data.1$PatientNum[annot_data.1$SegmentDisplayName == xl_annot$SegmentDisplayName]
#xl_annot$PatientNum <- as.factor(xl_annot$PatientNum)

colData <- xl_annot %>%
  dplyr::select(SegmentDisplayName, U6ATAC_RNAscopeClass, SegmentLabel) %>%
  distinct()

colData$condition <- with(colData, paste(U6ATAC_RNAscopeClass, SegmentLabel, sep = "_"))
rownames(colData) <- colData$SegmentDisplayName

# Prepare countData matrix
countData <- xl_counts
countData <- as.data.frame(countData)
rownames(countData) <- countData[, 1]
countData <- countData[, -1]

# Ensure the correct order of samples in countData
correct_order <- colData$SegmentDisplayName
countData <- countData[, correct_order]
stopifnot(all(colnames(countData) == correct_order))  # This should return TRUE

# Convert to matrix and round values
countData <- as.matrix(countData)
countData <- round(countData)

# Prepare DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ U6ATAC_RNAscopeClass)
dds <- DESeq(dds)  # Running DESeq2 analysis

