library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
#setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD")
processed_data <- as.data.frame(fread("TCGA_LUAD_ProcessedC_varc.csv"),)
processed_data_phen <- as.data.frame(fread("TCGA_LUAD_TPM_Regression.csv"))
processed_data <- data.frame(processed_data[,-1], row.names = processed_data[,1])
phen <- processed_data_phen$phen
values <- c()

for(i in 1:30){
sample <- sample(processed_data, 3000)
sample <- cbind(processed_data_phen$V1, sample, processed_data_phen$phen)
colnames(sample)[ncol(sample)] <- "phen"

normal_phen <- subset(sample, phen == "Solid Tissue Normal")
normal_phen <- subset(normal_phen, select = -phen)
normal_phen <- data.frame(normal_phen[,-1], row.names = normal_phen[,1])

affected_phen <- subset(sample, phen == "Primary Tumor")
affected_phen <- subset(affected_phen, select = -phen)
affected_phen <- data.frame(affected_phen[,-1], row.names = affected_phen[,1])

norm_corr <- cor(as.matrix(normal_phen), method = "pearson")
affected_corr <- cor(as.matrix(affected_phen), method = "pearson")

flattenCorrMatrix <- function(cormat){
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor = (cormat)[ut])
}    

norm_corr_flattened <- flattenCorrMatrix(norm_corr)
affected_corr_flattened <- flattenCorrMatrix(affected_corr)
corr_mat <- as.data.frame(cbind(norm_corr_flattened, affected_corr_flattened[,3]))
colnames(corr_mat) <- c("Gene 1", "Gene 2", "Normal Corr", "Affected Corr")
corr_mat$Diff <- abs(corr_mat$`Normal Corr` - corr_mat$`Affected Corr`)

df_ordered <- corr_mat[order(corr_mat$Diff),]
df_ordered <- df_ordered[!is.na(df_ordered$Diff),]
obs_sig_num <- floor(0.95 * length(df_ordered$Diff))
critical_value <- df_ordered[obs_sig_num, 5]
print(critical_value)
values <- c(values, critical_value)
}
print(mean(values))