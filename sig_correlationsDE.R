library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)

setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
sig_genes <- as.list(fread("Sig_DE_final.csv"))
processed_data <- as.data.frame(fread("TCGA_LUAD_TPM_Regression.csv"))
phen <- processed_data$phen
rownames(processed_data) <- processed_data[,1]
processed_data <- processed_data[,-1]
processed_data <- subset(processed_data, select = -phen)

index <- c()
all_genes <- colnames(processed_data)
for (i in 1:length(all_genes)){
  if(all_genes[i] %in% sig_genes$sig_genes){
    index <- c(index, i)
  }
}
sigDE <- subset(processed_data, select = index)
#Gene 1, Gene 2, Correlation
table_cor <- cor(sigDE, processed_data, method = "pearson")

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = colnames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}
flat_cor <- flattenCorrMatrix(table_cor)
#-0.12582, 0.59135
flat_cor <- flat_cor[flat_cor$cor <= -0.12582 | flat_cor$cor >= 0.59135,]

fwrite(as.data.frame(flat_cor), "SigCor_DE_final.csv", row.names = FALSE, col.names = TRUE)
