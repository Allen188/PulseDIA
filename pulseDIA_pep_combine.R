rm(list = ls())
req.pcg <- function(pcg){
  new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, dependencies = T)
  sapply(pcg, require, ch = T)
}
all.pcg <- c("stringr")
req.pcg(all.pcg)
sapply(all.pcg, library, character.only = T)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#----------------------------- Please modify the parameters ---------------------------------------------------

input_name <- "peptides_2019513.txt"
out_name <- "peptides_201cc9513_sum.txt"

#---------------------------------------------------------------------------------------------------------------

library(stringr)
df <- read.table(input_name,sep="\t",header = T,stringsAsFactors=F)
nm <- as.character(sapply(colnames(df)[-(1:2)],function(v) {str_sub(v,1,(regexpr("_part",v)-1))}))
data <- df
names(data)[3:ncol(df)] <- nm 
data1 <-  data.frame(peptide_group_label=data$peptide_group_label,prot=data$prot, data[,order(nm)+2],check.names = F)


df0 <- data1[,-c(1,2)]
df0 <- t(df0)
df0 <- data.frame(df0)
df0$label <- rownames(df0)
df0 <- df0[,c(ncol(df0),1:(ncol(df0)-1))]
df0$label <- sapply(df0$label,function(v){strsplit(v,"\\.")[[1]][1]})

result <- c()
k0 <- unique(df0$label)
for(lbl in k0){
  k1 <- df0[df0$label==lbl,]
  k2 <- apply(k1[,-1],2,function(v){mean(v,na.rm = T)})
  result <- rbind(result,k2)
}


rownames(result) <- k0

colnames(result) <- data$peptide_group_label
result <- data.frame(t(result))
clnames <- colnames(result)
result$peptide_group_label <- data$peptide_group_label
result$prot <- data$prot
result <- result[,c('peptide_group_label','prot',clnames)]
rownames(result) <- 1:dim(result)[1]
result1 <-  data.frame(peptide_group_label=data$peptide_group_label,prot=data$prot, result[,order(str_sub(names(result)[3:ncol(result)],-1,-1))+2])
##result1 <- result1[,-ncol(result1)]
write.table(result1,file=out_name,sep="\t",col.names = T,row.names = F,quote = F)