rm(list = ls())
req.pcg <- function(pcg){
  new <- pcg[!(pcg %in% installed.packages()[, "Package"])]
  if (length(new)) install.packages(new, dependencies = T)
  sapply(pcg, require, ch = T)
}
all.pcg <- c("readxl","stringr","magrittr")
req.pcg(all.pcg)
sapply(all.pcg, library, character.only = T)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#----------------------------- Please modify the parameters ---------------------------------------------------

wins_type <- "Fixed"                  #Fixed or Variable

pulse_num  <-   4
out_name <- "BRP_4_part"
out_filename_MS <- "4/MS" 
out_filename_tmp <- "4/Cirt"
overlap <-   FALSE                       # TRUE or FALSE

mz_start   <-   400
mz_end     <-   1200 
win_num    <-   100
input_name <- "BRP_QE_peptides.txt"


#----------------------------------------------------------------------------------------
df <- read.table(input_name,sep='\t',header=T,stringsAsFactors=F) %>% as.data.frame()

if(wins_type=="Variable" ){    
    outdata <- function(m,z,intensity){
      for (i in 1:length(z)) {
        if(stringr::str_detect(z[i],";")){
          str <- unlist(str_split(z[i], ";"))
          z[i] <- str[1]
          for (t in 2:length(str)) {
            z <- append(z,str[t])
            m <- append(m,m[i])
            intensity <- append(intensity,intensity[i])
          }
        }
      }
      mz <-  c(as.numeric(m)/as.numeric(z)) %>% round(1)
      data <- data.frame(mz=mz,intensity=intensity)
      return(data)
    }
    
    df1 <- outdata(df$Mass,df$Charges,df$Intensity)
    df1 <- na.omit(df1)
    
    df2 <- df1[order(df1$mz),]
    df2 <- df2[df2$mz>=mz_start & df2$mz<=mz_end,]
    row.names(df2) <- c(1:nrow(df2))
    df3 <- tapply(df2$intensity,factor(df2$mz),sum) %>% as.matrix()
    
    
    out_win <- function(){
      
        average <- sum(df3)/win_num/pulse_num
        for (bb in c(seq(0.01,0.3,0.01))) {
          group <- matrix(NA,ncol = 3,nrow = win_num*pulse_num)
          name <- row.names(df3)
          end=0
          for (gr in 1:(win_num*pulse_num)) {
            start=end+1
            if(start<=nrow(df3)){
              sum=0
              for (x in start:nrow(df3)) {
                sum=df3[x]+sum
                if(x==nrow(df3)){
                  end=x
                  break
                }
                if(sum >= average){
                  if(sum <= average*(1+bb)){
                    end=x
                    break
                  }
                  if(sum > average*(1+bb)){
                    if(sum==df3[x]){
                      end=x
                      sum=df3[x]
                      break
                    }else{
                      end=x-1
                      sum=sum-df3[x]
                      break
                    }
                  }
                }else{next}  
              }
              if(sum>0){
                group[gr,] <- c(name[start],name[end],sum)
              }
            }
          }
          if(name[end]==name[nrow(df3)]){
            break
          }
        }
          return(group)
    }
    
    outgr_no_ovlap <- function(){
      win <- out_win()
      lap <- matrix(NA,ncol = 3,nrow = win_num*pulse_num)
      for (i in 1:nrow(win)) {
        num1 <- as.numeric(win[i,1])
        num2 <- as.numeric(win[i,2])
        average <- num2-num1
          lap[i,] <- c(num1-0.5,num2+0.5,average+1)
          if(i==nrow(win))
          {
            lap[nrow(win),] <- c(num1-0.5,num2,average+0.5)
          }
      }
      return(lap)
    }
    
    outgr_ovlap <- function(){
      win <- out_win()
      lap <- matrix(NA,ncol = 3,nrow = win_num*pulse_num+1)
      for (i in 1:nrow(win)) {
        num1 <- as.numeric(win[i,1])
        num2 <- as.numeric(win[i,2])
        average <- num2-num1
          lap[i,] <- c(num1-average/2,num2+average/2,num2+average-num1)
          if(i==1 )
          {
            lap[1,] <- c(num1,num2+average/2,num2-num1+average/2)
          }
          if(i==nrow(win) )
          {
            lap[nrow(win),] <- c(num1-average/2,num2,num2-num1+average/2)
          }
      }
      lap[win_num*pulse_num+1,] <- c(lap[1,1],lap[1,1]+lap[1,3]/2,lap[1,3]/2)
      return(lap)
    }
    cc <- outgr_no_ovlap()
    sum(!is.na(cc[,2]))
    out_group <- function(){
      if(overlap==FALSE){ group <- outgr_no_ovlap()
      }else{
        group <- outgr_ovlap()}
      
      for (i in 1:sum(!is.na(group[,2]))) {
        if(group[i,3]<2){
          sum <- group[i,1]+group[i,2]
          group[i,1] <- round((sum)/2,1)-1
          group[i,2] <- round((sum)/2,1)+1
          group[i,3] <- 2
        }
      }
      title <- c("start","end","width")
      if(file.exists(out_filename_tmp)){
        unlink(out_filename_tmp, recursive=TRUE)
      } 
      dir.create(path=out_filename_tmp,recursive = TRUE)
      if(overlap==FALSE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          write.table(gp,paste0(out_filename_tmp,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title)
        }
      }
      if(overlap==TRUE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          if(j==1){
            front <- gp[nrow(gp),]
            gp <- gp[-nrow(gp),]
          }
          if(j==pulse_num){
            tm <- matrix(NA,ncol = 3,nrow = nrow(gp)+1)
            tm[2:(nrow(gp)+1),] <- gp
            tm[1,] <- front
            gp <- tm
          }
          write.table(gp,paste0(out_filename_tmp,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title)
        }
      }
    }
    
    out_dia <- function(){
      if(overlap==FALSE){ 
        group <- outgr_no_ovlap()
      }else{
        group <- outgr_ovlap()}
      
      for (i in 1:sum(!is.na(group[,2]))) {
        if(group[i,3]<2){
          sum <- group[i,1]+group[i,2]
          group[i,1] <- round((sum)/2,1)-1
          group[i,2] <- round((sum)/2,1)+1
          group[i,3] <- 2
        }
      }
      title <- c("Mass [m/z]","Formula [M]","CS [z]","Polarity","Start [min]","End [min]","NCE","Comment")
      if(file.exists(out_filename_MS)){
        unlink(out_filename_MS, recursive=TRUE)
      } 
      dir.create(path=out_filename_MS,recursive = TRUE)
      if(overlap==FALSE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          gp1 <- matrix(NA,ncol = 8,nrow = nrow(gp))
          gp1[,1] <- (gp[,1]+gp[,2])/2
          gp1[,4] <- c(rep("Positive",nrow(gp)))
          write.table(gp1,paste0(out_filename_MS,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title,na = "")
        }
      }
      if(overlap==TRUE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          if(j==1){
            front <- gp[nrow(gp),]
            gp <- gp[-nrow(gp),]
          }
          if(j==pulse_num){
            tm <- matrix(NA,ncol = 3,nrow = nrow(gp)+1)
            tm[2:(nrow(gp)+1),] <- gp
            tm[1,] <- front
            gp <- tm
          }
          gp1 <- matrix(NA,ncol = 8,nrow = nrow(gp))
          gp1[,1] <- (gp[,1]+gp[,2])/2
          gp1[,4] <- c(rep("Positive",nrow(gp)))
          write.table(gp1,paste0(out_filename_MS,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title,na = "")
        }
      }
    }
    out_group()
    out_dia()
}else
  if(wins_type=="Fixed" ){
    
    out_win <- function(){
      group <- matrix(NA,ncol = 2,nrow = 24)
      start=400
      for (gr in 1:20) {
        group[gr,] <- c(start,start+20)
        start=start+20
      }
      group[21,] <- c(800,860)
      group[22,] <- c(860,940)
      group[23,] <- c(940,1060)
      group[24,] <- c(1060,1200)
      return(group)
    }
    xx <- out_win()
    outgr_no_ovlap <- function(){
      win <- out_win()
      lap <- matrix(NA,ncol = 3,nrow = 24*pulse_num)
      for (i in 1:nrow(win)) {
        num1 <- as.numeric(win[i,1])
        num2 <- as.numeric(win[i,2])
        average <- ceiling((num2-num1)/pulse_num *10)/10
        for (e in 1:pulse_num) {
          lap[pulse_num*(i-1)+e,] <- c(num1-0.5,num1+average+0.5,average+1)
          if(i==nrow(win) & e==pulse_num)
          {
            lap[pulse_num*(i-1)+e,] <- c(num1-0.5,num1+average,average+0.5)
          }
          num1 <- num1+average
        }
      }
      return(lap)
    }
    
    outgr_ovlap <- function(){
      win <- out_win()
      lap <- matrix(NA,ncol = 3,nrow = 24*pulse_num+1)
      for (i in 1:nrow(win)) {
        num1 <- as.numeric(win[i,1])
        num2 <- as.numeric(win[i,2])
        average <- ceiling((num2-num1)/pulse_num *10)/10
        for (e in 1:pulse_num) {
          lap[pulse_num*(i-1)+e,] <- c(num1,num1+2*average,2*average)
          if(i==nrow(win) & e==pulse_num)
          {
            lap[pulse_num*(i-1)+e,] <- c(num1,num1+average,average)
          }
          num1 <- num1+average
        }
      }
      lap[24*pulse_num+1,] <- c(lap[1,1],lap[1,1]+lap[1,3]/2,lap[1,3]/2)
      return(lap)
    }
    cc <- outgr_ovlap()
    out_group <- function(){
      if(overlap==FALSE){ group <- outgr_no_ovlap()
      }else{
        group <- outgr_ovlap()}
      
      for (i in 1:nrow(group)) {
        if(group[i,3]<2){
          group[i,1] <- group[i,1]-1
          group[i,2] <- group[i,2]+1
          group[i,3] <- group[i,3]+2
        }
      }
      title <- c("start","end","width")
      if(file.exists(out_filename_tmp)){
        unlink(out_filename_tmp, recursive=TRUE)
      } 
      dir.create(path=out_filename_tmp,recursive = TRUE)
      if(overlap==FALSE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          write.table(gp,paste0(out_filename_tmp,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title)
        }
      }
      if(overlap==TRUE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          if(j==1){
            front <- gp[nrow(gp),]
            gp <- gp[-nrow(gp),]
          }
          if(j==pulse_num){
            tm <- matrix(NA,ncol = 3,nrow = nrow(gp)+1)
            tm[2:(nrow(gp)+1),] <- gp
            tm[1,] <- front
            gp <- tm
          }
          write.table(gp,paste0(out_filename_tmp,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title)
        }
      }
    }
    
    out_dia <- function(){
      if(overlap==FALSE){ 
        group <- outgr_no_ovlap()
      }else{
        group <- outgr_ovlap()}
      
      for (i in 1:nrow(group)) {
        if(group[i,3]<2){
          group[i,1] <- group[i,1]-1
          group[i,2] <- group[i,2]+1
          group[i,3] <- group[i,3]+2
        }
      }
      title <- c("Mass [m/z]","Formula [M]","CS [z]","Polarity","Start [min]","End [min]","NCE","Comment")
      if(file.exists(out_filename_MS)){
        unlink(out_filename_MS, recursive=TRUE)
      } 
      dir.create(path=out_filename_MS,recursive = TRUE)
      if(overlap==FALSE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          gp1 <- matrix(NA,ncol = 8,nrow = nrow(gp))
          gp1[,1] <- (gp[,1]+gp[,2])/2
          gp1[,4] <- c(rep("Polarity",nrow(gp)))
          write.table(gp1,paste0(out_filename_MS,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title,na = "")
        }
      }
      if(overlap==TRUE){
        for (j in 1:pulse_num) {
          gp <- group[seq(j,nrow(group),pulse_num),]
          if(j==1){
            front <- gp[nrow(gp),]
            gp <- gp[-nrow(gp),]
          }
          if(j==pulse_num){
            tm <- matrix(NA,ncol = 3,nrow = nrow(gp)+1)
            tm[2:(nrow(gp)+1),] <- gp
            tm[1,] <- front
            gp <- tm
          }
          gp1 <- matrix(NA,ncol = 8,nrow = nrow(gp))
          gp1[,1] <- (gp[,1]+gp[,2])/2
          gp1[,4] <- c(rep("Polarity",nrow(gp)))
          write.table(gp1,paste0(out_filename_MS,"/",out_name,"_",j,".csv"), row.names = FALSE,sep = ",", col.names = title,na = "")
        }
      }
    }
    out_group()
    out_dia()
  }else{
    print("Error: win_type is wrong!!!")
  }







