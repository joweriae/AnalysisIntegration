library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
#setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD")
processed_data <- as.data.frame(fread("TCGA_LUAD_ProcessedC_varc.csv"))
processed_data <- data.frame(processed_data[,-1], row.names = processed_data[,1])
processed_data_phen <- as.data.frame(fread("TCGA_LUAD_TPM_Regression.csv"))
sample_data <- as.data.frame(fread("clinical.cart.2020-10-15/clinical.tsv"))
phen <- processed_data_phen$phen
sample_names <- processed_data[1]

sex_table <- matrix(data = NA, nrow = 0, ncol = 2)
colnames(sex_table) <- c("Sample ID", "Sex")
for (i in 1:length(sample_names[,1])){
  name <- substr(sample_names[i, 1], 1, 12)
  index <- match(name, sample_data$case_submitter_id)
  sex <- sample_data[index, 12]
  sex_table <- rbind(sex_table, c(name, sex))
}
processed_data$sex <- sex_table[,2]
females <- subset(processed_data, sex == "female")
males <- subset(processed_data, sex == "male")

normal_phen <- subset(females, phen == "Solid Tissue Normal")
normal_phen <- subset(normal_phen, select = -phen)
normal_phen <- subset(normal_phen, select = -sex)
normal_phen <- data.frame(normal_phen[,-1], row.names = normal_phen[,1])


affected_phen <- subset(females, phen == "Primary Tumor")
affected_phen <- subset(affected_phen, select = -phen)
affected_phen <- subset(affected_phen, select = -sex)
affected_phen <- data.frame(affected_phen[,-1], row.names = affected_phen[,1])

flattenCorrMatrix <- function(cormat) {
  ut <- upper.tri(cormat)
  data.table(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  = (cormat)[ut]
  )
}
norm_corr <- cor(as.matrix(normal_phen), method = "pearson")
affected_corr <- cor(as.matrix(affected_phen), method = "pearson")

norm_corrs <- flattenCorrMatrix(norm_corr)
aff_corrs <- flattenCorrMatrix(affected_corr)

merged <- merge(norm_corrs, aff_corrs, by.x=c("row", "column"), by.y=c("row", "column"))
merged$diff <- abs(merged$cor.x - merged$cor.y)
sig_cor <- merged[merged$diff >= 0.4981441,]

fwrite(as.data.frame(sig_cor), "Sig_Correlations_females.csv", row.names = FALSE, col.names = TRUE)

normal_phen <- subset(males, phen == "Solid Tissue Normal")
normal_phen <- subset(normal_phen, select = -phen)
normal_phen <- subset(normal_phen, select = -sex)
normal_phen <- data.frame(normal_phen[,-1], row.names = normal_phen[,1])


affected_phen <- subset(males, phen == "Primary Tumor")
affected_phen <- subset(affected_phen, select = -phen)
affected_phen <- subset(affected_phen, select = -sex)
affected_phen <- data.frame(affected_phen[,-1], row.names = affected_phen[,1])

norm_corr <- cor(as.matrix(normal_phen), method = "pearson")
affected_corr <- cor(as.matrix(affected_phen), method = "pearson")

norm_corrs <- flattenCorrMatrix(norm_corr)
aff_corrs <- flattenCorrMatrix(affected_corr)

merged <- merge(norm_corrs, aff_corrs, by.x=c("row", "column"), by.y=c("row", "column"))
merged$diff <- abs(merged$cor.x - merged$cor.y)
sig_cor <- merged[merged$diff >= 0.4981441,]

fwrite(as.data.frame(sig_cor), "Sig_Correlations_males.csv", row.names = FALSE, col.names = TRUE)
