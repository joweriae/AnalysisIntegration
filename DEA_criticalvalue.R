library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
#setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD")

processed_data_phen <- as.data.frame(fread("TCGA_LUAD_TPM_Regression.csv"))
phen <- processed_data_phen$phen
processed_data_phen <- subset(processed_data_phen, select = -phen)

values <- c()
for(i in 1:10500){
  sample <- sample(processed_data_phen, 2)
  #sample <- cbind(processed_data_phen$V1, sample, processed_data_phen$phen)
  #colnames(sample)[ncol(sample)] <- "phen"
  
  rand_corr <- cor(sample[,1], sample[,2], method = "pearson")
  values <- c(values, rand_corr)
}
values <- values[complete.cases(values)]
val_ordered <- values[order(values)]
obs_sig_num <- length(values) - floor(0.975 * length(values))
sig_num_pos <- val_ordered[length(val_ordered)-obs_sig_num] #0.5913525
sig_num_neg <- val_ordered[obs_sig_num] #-0.1258232