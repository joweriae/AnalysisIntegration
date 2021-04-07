library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
#setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD")
processed_data <- as.data.frame(fread("TestDataset.csv"))

normal_phen <- subset(processed_data, phen == "Solid Tissue Normal")
normal_phen <- subset(normal_phen, select = -phen)
normal_phen_rownames <- data.frame(normal_phen[,-1], row.names = normal_phen[,1])
normal_phen_transposed <- transpose(normal_phen_rownames)
rownames(normal_phen_transposed)<-colnames(normal_phen_rownames)
colnames(normal_phen_transposed)<-rownames(normal_phen_rownames)
normal_phen_means <- rowMeans(normal_phen_transposed)

affected_phen <- subset(processed_data, phen == "Primary Tumor")
affected_phen <- subset(affected_phen, select = -phen)
affected_phen_rownames <- data.frame(affected_phen[,-1], row.names = affected_phen[,1])
affected_phen_transposed <- transpose(affected_phen_rownames)
rownames(affected_phen_transposed)<-colnames(affected_phen_rownames)
colnames(affected_phen_transposed)<-rownames(affected_phen_rownames)
affected_phen_means <- rowMeans(affected_phen_transposed)
comparitive_expr <- cbind(normal_phen_means, affected_phen_means)

#geneIDs <- row.names(normal_phen_transposed)
#cor_mat_norm <- matrix(, nrow = length(geneIDs), ncol = length(geneIDs))
#p_mat_norm <- matrix(, nrow = length(geneIDs), ncol = length(geneIDs))
#rownames(cor_mat_norm) <- geneIDs
#colnames(cor_mat_norm) <- geneIDs
#for(id in length(geneIDs)){
#  for(id2 in length(geneIDs)){
#    cor_test <- cor.test(as.numeric(data.matrix(normal_phen_rownames[id])), as.numeric(data.matrix(normal_phen_rownames[id2])), method = "pearson")
#    cor_mat_norm[id, id2] <- cor_test$estimate
#    p_mat_norm[id, id2] <- cor_test$p.value

#  }
#}
#rownames(cor_mat_norm) <- geneIDs
#colnames(cor_mat_norm) <- geneIDs
#rownames(p_mat_norm) <- geneIDs
#colnames(p_mat_norm) <- geneIDs

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
# flattenCorrMatrix <- function(cormat, pmat) {
#   ut <- upper.tri(cormat)
#   data.frame(
#     row = rownames(cormat)[row(cormat)[ut]],
#     column = rownames(cormat)[col(cormat)[ut]],
#     cor  =(cormat)[ut],
#     p = pmat[ut]
#   )
# }

norm_corr <- cor(as.matrix(normal_phen_rownames), method = "pearson")
affected_corr <- cor(as.matrix(affected_phen_rownames), method = "pearson")
#formatted <- flattenCorrMatrix(cor_mat_norm$r, corr_mat_norm$P)

p_val_calc <- function(r){
  z <- 0.5 * log((1+r)/(1-r))
  zse <- 1/sqrt(54-3)
  min(pnorm(z, sd = zse), pnorm(z, lower.tail = F, sd = zse))*2
}

p_val_data_norm <- apply(norm_corr, 2, Vectorize(p_val_calc))
p_val_data_aff <- apply(affected_corr, 2, Vectorize(p_val_calc))

sig_norm <- which(p_val_data_norm <= 0.05, arr.ind = T)
sig_aff <- which(p_val_data_aff <= 0.05, arr.ind = T)
names <- row.names(norm_corr)
names_aff <- row.names(affected_corr)
cor_mat_norm <- matrix(data = NA, nrow = nrow(sig_norm), ncol = ncol(sig_norm))
cor_mat_aff <- matrix(data = NA, nrow = nrow(sig_aff), ncol = ncol(sig_aff))
for(row in 1:nrow(sig_norm)){
  for(col in 1:ncol(sig_norm)){
    cor_mat_norm[row, col] <- names[sig_norm[row, col]]
  }
}
for(row in 1:nrow(sig_aff)){
  for(col in 1:ncol(sig_aff)){
    cor_mat_aff[row, col] <- names_aff[sig_aff[row, col]]
  }
}

fwrite(as.data.frame(cor_mat_norm), "TestOutput.csv", row.names = TRUE, col.names = TRUE)
