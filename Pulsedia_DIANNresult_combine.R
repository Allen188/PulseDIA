## PulseDIA peptide result combine
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

file <- list.files("./","pep")
n=1
list <- list()
for (i in file) {
  list[[n]] <- read.table(i,sep = "\t",header = T,stringsAsFactors = F)
  n=n+1
}
tmp <- merge(list[[1]],list[[2]],by.y = c("pep","prot"),all=T)
if(length(file)>2){
  for (ii in 3:length(file)) {
    tmp <- merge(tmp,list[[ii]],by.y = c("pep","prot"),all=T)
  }
}

overlap_calcu <- "mean"       ## "mean"  or  "max"
out_name <- "diann_cca_pep20191206_combine.txt"

#---------------------------------------------------------------------------------------------------------------

library(stringr)
df <- tmp
nm <- as.character(sapply(colnames(df)[-(1:2)],function(v) {str_split(v,"\\.")[[1]][6]}))
nm <- as.character(sapply(nm,function(v){str_split(v,"_part")[[1]][1]}))
data <- df
names(data)[3:ncol(df)] <- nm 
data1 <-  data.frame(peptide_group_label=data[,1],prot=data[,2], data[,order(nm)+2],check.names = F)


df0 <- data1[,-c(1,2)]
df0 <- t(df0)
df0 <- data.frame(df0)
df0$label <- rownames(df0)
df0 <- df0[,c(ncol(df0),1:(ncol(df0)-1))]
df0$label <- sapply(df0$label,function(v){strsplit(v,"\\.")[[1]][1]})

result <- c()
k0 <- unique(df0$label)
if(overlap_calcu=="mean"){
  for(lbl in k0){
    k1 <- df0[df0$label==lbl,]
    k2 <- apply(k1[,-1],2,function(v){mean(v,na.rm = T)})
    result <- rbind(result,k2)
  }
}

if(overlap_calcu=="max"){
  for(lbl in k0){
    k1 <- df0[df0$label==lbl,]
    k2 <- apply(k1[,-1],2,function(v){ifelse(sum(!is.na(v))>0,max(v,na.rm = T),NA) })
    result <- rbind(result,k2)
  }
}

rownames(result) <- k0

colnames(result) <- data1$peptide_group_label
result <- data.frame(t(result))
clnames <- colnames(result)
result$peptide_group_label <- data1$peptide_group_label
result$prot <- data1$prot
result <- result[,c('peptide_group_label','prot',clnames)]
rownames(result) <- 1:dim(result)[1]
result1 <-  data.frame(peptide_group_label=data1$peptide_group_label,prot=data1$prot, result[,order(str_sub(names(result)[3:ncol(result)],-1,-1))+2])
##result1 <- result1[,-ncol(result1)]
write.table(result1,file=out_name,sep="\t",col.names = T,row.names = F,quote = F)


#########################################################################################
## PulseDIA protein result combine

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
file <- list.files("./","prot")
n=1
list <- list()
for (i in file) {
  list[[n]] <- read.table(i,sep = "\t",header = T,stringsAsFactors = F)
  n=n+1
}
tmp <- merge(list[[1]],list[[2]],by.y = c("prot"),all=T)
if(length(file)>2){
  for (ii in 3:length(file)) {
    tmp <- merge(tmp,list[[ii]],by.y = c("prot"),all=T)
  }
}

overlap_calcu <- "mean"       ## "mean"  or  "max"
out_name <- "diann_cca_prot20191206_combine.txt"

#---------------------------------------------------------------------------------------------------------------

library(stringr)
df <- tmp
nm <- as.character(sapply(colnames(df)[-1],function(v) {str_split(v,"\\.")[[1]][6]}))
nm <- as.character(sapply(nm,function(v){str_split(v,"_part")[[1]][1]}))
data <- df
names(data)[2:ncol(df)] <- nm 
data1 <-  data.frame(prot=data[,1], data[,order(nm)+1],check.names = F)


df0 <- data1[,-1]
df0 <- t(df0)
df0 <- data.frame(df0)
df0$label <- rownames(df0)
df0 <- df0[,c(ncol(df0),1:(ncol(df0)-1))]
df0$label <- sapply(df0$label,function(v){strsplit(v,"\\.")[[1]][1]})

result <- c()
k0 <- unique(df0$label)
if(overlap_calcu=="mean"){
  for(lbl in k0){
    k1 <- df0[df0$label==lbl,]
    k2 <- apply(k1[,-1],2,function(v){mean(v,na.rm = T)})
    result <- rbind(result,k2)
  }
}

if(overlap_calcu=="max"){
  for(lbl in k0){
    k1 <- df0[df0$label==lbl,]
    k2 <- apply(k1[,-1],2,function(v){ifelse(sum(!is.na(v))>0,max(v,na.rm = T),NA) })
    result <- rbind(result,k2)
  }
}

rownames(result) <- k0

colnames(result) <- data1$prot
result <- data.frame(t(result))
clnames <- colnames(result)
result$prot <- data1$prot
result <- result[,c('prot',clnames)]
rownames(result) <- 1:dim(result)[1]
result1 <-  data.frame(prot=data1$prot, result[,order(str_sub(names(result)[2:ncol(result)],-1,-1))+1])
##result1 <- result1[,-ncol(result1)]
write.table(result1,file=out_name,sep="\t",col.names = T,row.names = F,quote = F)
