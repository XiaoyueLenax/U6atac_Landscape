setwd("/Users/lucaroma/Desktop/smarca2_4_revision/msk_ar_signature/")
library(ggplot2)
library(ggsignif)
vector <- read.table("LnCaP.reference.AR_score.txt", sep = "\t",header = T)

#input <- read.table("input_table.txt",sep = "\t", header = T)
#msk_samples_log_nepc

y <- read.table("AR_list.txt",header = T,sep = "\t")

y=DGEList(counts=y[,-1],genes=y[,1])



logFPKM <- log2(y$counts+1)
#logFPKM <- as.data.frame(logFPKM)

z <- ""
for (i in colnames(logFPKM)){
  a <- cor.test(logFPKM[,i],vector$value, method = "pearson") # check that the logFPKM has the same gene order of reference vector
  z[i] = a$estimate
}