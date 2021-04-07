library(data.table)
library(stringr)
setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')

#TPM Regression doc with phenotypic data
processed_data <- as.data.frame(fread("TCGA_LUAD_TPM_Regression.csv"))
processed_data <- data.frame(processed_data[,-1], row.names = processed_data[,1])
#Save the phenotypic data seperately
phen <- processed_data$phen
#Remove the phenotype data from the dataset we will be using
processed_data <- as.matrix(subset(processed_data, select = -phen))
#Calculate the variance for all
variances <- apply(processed_data, 2, var)

#Empty variance matrix
var_probe <- matrix(data = NA, nrow = 0, ncol = 2)
var_probe <- cbind(colnames(processed_data), variances)
colnames(var_probe) <- c("Gene Name", "Variance")

var_ordered <- as.data.frame(var_probe)
var_ordered <- var_ordered[order(var_ordered$Variance),]
var_ordered <- subset(var_ordered, var_ordered$Variance >= 1)
#num <- nrow(var_ordered)
#cutoff <- ceiling(num * 0.25)
#rem <- 1:cutoff
#remove <- var_ordered[rem,][,1]
keep <- var_ordered[,1]
complete <- subset(processed_data, select = keep)
fwrite(as.data.table(keep), "kept_genes.csv")
complete <- cbind(rownames(complete), complete, phen)
fwrite(as.data.table(complete), "TCGA_LUAD_ProcessedC_varc.csv", row.names = TRUE)
