library(data.table)
library(stringr)
library(Hmisc)
library(rstatix)
library(WebGestaltR)
#setwd('~//..//..//work//long_lab//joweria//TCGA-LUAD')
setwd("C://Users//babyb//OneDrive//Desktop//TCGA-LUAD//DEA March 7 2021")
sig_cor <- as.data.frame(fread("SigCor_DE_final.csv"))
sig_genes <- as.data.frame(fread("Sig_DE_final.csv"))
colnames(sig_cor) <- c("Gene 1", "Gene 2", "Sig Cor")
cors <- c(sig_cor[,1], sig_cor[,2])
#count_cors <- aggregate(data.frame(count = cors), list(value = cors), length)
#ordered_cors <- count_cors[order(count_cors$count, decreasing = TRUE),]
count <- aggregate(sig_cor[,1], list(value = sig_cor[,1]), length)
ord <- count[order(count[,2], decreasing = TRUE),]
fwrite(as.data.frame(ord), "Cor_count_DE_final.csv", row.names = FALSE, col.names = TRUE)
ord <- as.data.frame(fread("Cor_count_DE_final.csv"))
thres_cors <- ord[ord$x >= 500,]
these <- ord$value
fin_sig <- as.data.frame(fread("Cor_count_DE_final.csv"))
#thres_cors <- as.matrix(fin_sig)
class(thres_cors[,2]) <- "numeric"
for(val in 1:length(thres_cors[,2])){
  value <- thres_cors[val, 1]
  sig <- sig_cor[sig_cor$`Gene 1` == value,]
  sig2 <- sig_cor[sig_cor$`Gene 2` == value,]
  sig3 <- test3 <-  rbind(sig, sig2)
  print(sig3[4,])
}
pathways <- matrix(data = NA, ncol = 2)
for (i in 1:length(thres_cors$x)){
  print(i)
  print(length(thres_cors$x))
  x <- thres_cors$value[i]
  test <- sig_cor[sig_cor$`Gene 1` == x,]
  test_proc <- test[,2]
  test2 <- sig_cor[sig_cor$`Gene 2` == x,]
  test2_proc <- test2[,1]
  c1 <- c(test_proc, test2_proc)
  c1 <- c1[c1 != x]
  string <- paste(x, "Enrich.txt", sep = "")
fwrite(as.list(c1), string, row.names = FALSE, col.names = FALSE, sep = "\n")
  enrichment_results <- WebGestaltR(interestGene = c1, organism = "hsapiens", enrichDatabase="pathway_KEGG", 
                                    interestGeneType="ensembl_gene_id", referenceSet = "genome",
                                    referenceGeneType = "ensembl_gene_id", isOutput = TRUE, projectName = string)
  len <- length(enrichment_results$description)
  for (j in 1:len){
    desc <- enrichment_results$description[j]
    if(!length(desc)){
      next
    }
    if(desc %in% pathways){
      position <- which(pathways[,1] == desc)
      count <- as.numeric(pathways[position, 2])
      pathways[position, 2] <- count + 1
    }
    else{
      pathways <- rbind(pathways, c(desc, 1))
    }
  }
}
path <- pathways
path[,2] <- as.numeric(path[,2])

fwrite(as.data.frame(path), "pathways_DECor_final.csv", row.names = FALSE, col.names = FALSE)