rm(list = ls())
req.pcg <- function(pcg){
  new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, dependencies = T)
  sapply(pcg, require, ch = T)
}
all.pcg <- c("stringr","readxl")
req.pcg(all.pcg)
sapply(all.pcg, library, character.only = T)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#----------------------------- Please modify the parameters ---------------------------------------------------

in1 <- read.table("20191025_003142_F20191024qianlj_CRDD_WV_MLhighpHlib_part1.xls",sep = "\t",header = T,stringsAsFactors = F)
in2 <- read.table("20191025_091505_F20191024qianlj_CRDD_WV_MLhighpH_part2.xls",sep = "\t",header = T,stringsAsFactors = F)
# in3 <- read.table("20191021_194458_B20191016qianlj_m1A_TSLlib_part3.xls",sep = "\t",header = T,stringsAsFactors = F)
# in4 <- read.table("20191021_194511_B20191016qianlj_m1A_TSLlib_part4.xls",sep = "\t",header = T,stringsAsFactors = F)

out_name <- "20191025_091505_F20191024qianlj_CRDD_WV_MLhighpH_pep.txt"

#---------------------------------------------------------------------------------------------------------------
tmp <- merge(in1,in2,by.y = c("PG.ProteinNames","PEP.StrippedSequence"),all=T)


tmp1 <- tmp
tmp1[tmp1=="NaN"] <- NA
tmp1[tmp1=="Filtered"] <- NA
tmp1 <- as.data.frame(lapply(tmp1,as.numeric))

library(reshape2) 
tmp3 <- data.frame(pep=tmp$PEP.StrippedSequence,prot=tmp$PG.ProteinNames,tmp1[,3:ncol(tmp1)])


df <- tmp3
names(df)[1:2] <- c("peptide_group_label","prot")
nm <- as.character(sapply(colnames(df)[-(1:2)],function(v) {str_split(v,"_part")[[1]][1]}))
nm1 <- as.character(sapply(nm,function(v) {str_split(v,"CRDD_")[[1]][2]}))

data <- df
names(data)[3:ncol(df)] <- nm1
data1 <-  data.frame(peptide_group_label=data$peptide_group_label,prot=data$prot, data[,order(nm1)+2],check.names = F)


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

result1 <-  data.frame(peptide_group_label=data$peptide_group_label,prot=data$prot, result[,order(names(result)[3:ncol(result)])+2])

write.table(result1,file=out_name,sep="\t",col.names = T,row.names = F,quote = F)
