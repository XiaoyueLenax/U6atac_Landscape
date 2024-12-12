library(DESeq2)
library(tidyverse)
#BiocManager::install("dorothea")
#BiocManager::install("decoupleR")
BiocManager::install("viper")
library(dorothea)
library(viper)
library(decoupleR)
library(ggplot2)
library(dplyr)
library(viper)
# ========================= from site
net <- decoupleR::get_dorothea(levels = c('A', 'B', 'C', 'D'))
head(net)
n_genes <- net %>%
  group_by(source) %>%
  summarize(n = n())

ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")

n_edges <- net %>%
  group_by(confidence) %>%
  summarize(n = n())

ggplot(data=n_edges, aes(x=confidence, y=log10(n), color=confidence, fill=confidence)) +
  geom_bar(stat="identity") +
  theme(text = element_text(size=12)) +
  xlab('Confidence') +
  ylab('log10(Number of edges)') +
  theme_bw() +
  theme(legend.position = "none")

# ====================== explore TRI ==================================
net <- decoupleR::get_collectri(split_complexes = FALSE)
head(net)
n_genes <- net %>%
  group_by(source) %>%
  summarize(n = n())

ggplot(data=n_genes, aes(x=n)) +
  geom_density() +
  theme(text = element_text(size=12)) +
  xlab('Number of target genes') +
  ylab('densities') +
  theme_bw() +
  theme(legend.position = "none")
# =============================== dataset ======================================
e2f_list <- read.csv("dorothea_abc_tf_list.csv")
# Question: how to classify them into E2F targets? literature research?
# extract the TF that contains the E2F transcription factors. 
# Dorothea list is not on msigb, reformatted to be used as a ** GSCA ** 
# Can either run the whole dorothea pipeline or use viper to analyse 
# https://guolab.wchscu.cn/GSCA/#/
# https://saezlab.github.io/dorothea/reference/run_viper.html 
# gene mapping loss - check with online how can I solve 
# - as long as we dont lose the gene and are able to map the canonical ones, should be fine.
e2f_list
# manually check the failed mapped ones, and see if its a typo or I can manually change that
# in DESEQ2 - model patient variability to correct for batch effect among patients
# go back to the clinical data and check whether there is some variability among patients

# from pca: look at the genes that are contributing to PC1 and PC2. 
# 16 vs rest of them to explore what explains the difference. 
# speaks in terms on patient heterogeneity, patients cluster together no big
# variation within the patients themselves
# log transforms the visualization within PCA?
# proliferation score
# slides: as many high and as many low. less TME in the brain since sample was poorer.
# during selection, comparision between patient was not done. 
# match the U6atac images. 

# ----------------- with viper -----------------------------------
# Extract normalized counts
normalized_counts <- counts(dds, normalized=TRUE)