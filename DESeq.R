library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
library(DESeq2)
library(ggplot2)
library(corrr)
library("BiocParallel")


#setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD")
read_counts <- as.data.frame(fread("Count_Data.csv"))
remove_genes <- as.data.frame(fread("removed_genes.csv"))
kept_genes <- as.data.frame(fread("kept_genes.csv"))
colnames(read_counts)[1] <- "ensgene"

sample_info <- fread(list.files(pattern = "gdc_sample_sheet*", full.names = TRUE))
sample_info$`File Name` <- str_replace_all(sample_info$`File Name`, "\\.gz", "")
sample_info <- sample_info[endsWith(sample_info$`File Name`, "htseq.counts"),]
sample_info$`File Name` <- str_replace_all(sample_info$`File Name`, "\\.counts", "")

n <- colnames(read_counts)
#r <- t(read_counts)
#colnames(r)<-r[1,]
#r <- as.data.frame(r, row.names = TRUE)
#new1 <- matrix(nrow = 540, ncol = 0)
#for(x in 1:nrow(kept_genes)){
#  kept <- kept_genes[x, 1]
#  new1 <- cbind(new1, r[kept])
#}
#sample_name <- colnames(read_counts)
#rownames(new1) <- sample_name
#read_counts <- t(new1)
gene_names <- read_counts[,1]
phen <- c()
for(i in 2:length(n)){
  x <- match(n[i], sample_info$`Sample ID`)
  phen <- c(phen, sample_info$`Sample Type`[x])
}
meta <- cbind(n[-1], phen)
colnames(meta)[1]<- "id"

meta2 <- meta
meta_df <- as.data.frame(meta)
meta_df <- subset(meta_df, phen != "Recurrent Tumor")
read_counts <- subset(read_counts, select = -`TCGA-50-5946-02A`)
read_counts <- subset(read_counts, select = -`TCGA-50-5066-02A`)
rownames(read_counts) <- read_counts$ensgene
read_counts <- subset(read_counts, select = -`ensgene`)
read_counts_m <- as.matrix(read_counts)
class(read_counts_m) <- "numeric"
dds_rn <- DESeqDataSetFromMatrix(countData = read_counts_m, colData = meta_df, design = ~phen)
#
#keep <- rowSums(counts(dds)) > 1
#dds <- dds[keep,]
keep <- rowSums(counts(dds_rn)) >= 50
dds_rn <- dds_rn[keep,]

dds_rn <- DESeq(dds_rn, parallel = TRUE, BPPARAM = SnowParam(4))

resLFC1 <- results(dds_rn, lfcThreshold=2, alpha = 0.05)
res_sig <- resLFC1[complete.cases(resLFC1),]
res_sig <- res_sig[res_sig$padj < 0.01,]
sig_genes <- rownames(res_sig)
fwrite(as.data.frame(sig_genes), "Sig_DE_final.csv", row.names = FALSE, col.names = TRUE)

