library(data.table)
library(stringr)
library(rstatix)

setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
correlation_data <- as.matrix.data.frame(fread("TCGA_LUAD_cor.csv"))

p_val_calc <- function(r){
  z <- 0.5 * log((1+r)/(1-r))
  zse <- 1/sqrt(54-3)
  min(pnorm(z, sd = zse), pnorm(z, lower.tail = F, sd = zse))*2
}

correlation_data <- matrix(data = c(0.1, 0.2, 0.3), nrow = 3, ncol = 3)

p_val_data <- apply(correlation_data, 2, Vectorize(p_val_calc))

Which.names <- function(DF, value){
  ind <- which(DF<value, arr.ind=TRUE)
  paste(rownames(DF)[ind[1:nrow(ind)]],  colnames(DF)[ind[2]], sep=', ')
}

