# PulseDIA
A novel gas phase separation based MS acquisition method could improve the specificity and sensitivity of Orbitrap-based DIA analysis.
<br>
## Instruction
### PulseDIA_calcu_win.R
This script can generate a customized PulseDIA isolation window scheme.

#### Please set the required parameters according to experiments
```
#----------------------------- Please modify the parameters ---------------------------------------------------

wins_type <- "Fixed"                  #Fixed or Variable

pulse_num  <-   4
out_name <- "BRP_4_part"
out_filename_MS <- "4/MS" 
out_filename_Cirt <- "4/Cirt"
overlap <-   FALSE                       # TRUE or FALSE

mz_start   <-   400
mz_end     <-   1200 
win_num    <-   24
input_name <- "BRP_QE_peptides.txt"

#----------------------------------------------------------------------------------------------------------------
```
<br>wins_type: type of the designed window width (default: `Fixed` ; valid: `Fixed`, `Variable`), `Fixed` means the experiment with a fixed width window, `Variable` means the experiment with variable width windows that designed according to the ion density.
<br>pulse_num: number of injections
<br>out_name: the name of output file
<br>out_filename_MS: the output folder name. The generated files in this folder will be used for mass spectrometry
<br>out_filename_Cirt: the output folder name. The generated files in this folder will be used to build the windows file of CiRT for OpenSWATH analysis.
<br>overlap: Whether the windows have overlap(default: `FALSE` ;valid: `TRUE`, `FALSE`), `FALSE` means 1 Thomas overlap between two adjacent windows; `TRUE` means half width overlap between two adjacent windows
#### If you choose wins_type parameter to be `Variable`, you need to set the following parameters
<br>mz_start: The start m/z of MS1 acqusition range
<br>mz_end: The end m/z of MS1 acqusition range
<br>win_num: The number of isolation windows, default is 24 Windows
<br>input_name: precursor ion intensity file for special samples, window width designed according to the precursor ion density from this file

### PulseDIA_DIANNreport_extract.Rmd
This script is used to extract peptide or protein quantification results from DIA-NN output
```
df <- read.delim("report.tsv")
bb <- which(df$Precursor.Quantity==0)
df <- df[-bb,]
df1 <- data.frame(sample=df$File.Name,pep=df$Stripped.Sequence,prot=df$Protein.Ids,intesity=df$Precursor.Quantity)
df2 <- dcast(df1,pep+prot~sample,value.var = 'intesity',mean)
sum_dia <- apply(df2,2,function(x) sum(!is.na(x)))
write.table(df2,file="diann_CCA_pd2_pep20191206.txt",sep="\t",col.names = T,row.names = F,quote = F)
```
The code above is used to extract peptide quantitative results and export them
```
df <- read.delim("report.tsv")
bb <- which(df$Precursor.Quantity==0)
df <- df[-bb,]
#df <- df[-which(df$Protein.Ids==""),]
df1 <- data.frame(sample=df$File.Name,prot=df$Protein.Ids,intesity=df$PG.Quantity)
df1 <- unique(df1)
library(reshape2)
df2 <- dcast(df1,prot~sample,value.var = 'intesity',mean)
sum_dia <- apply(df2,2,function(x) sum(!is.na(x)))
write.table(df2,file="diann_CCA_pd2_prot20191206.txt",sep="\t",col.names = T,row.names = F,quote = F)
```
The code above is used to extract protein quantitative results and export them

### Pulsedia_DIANNresult_combine.R
The script could to combine peptides or proteins quantitative results from multiple PulseDIA injections for the same sample.
#### Please set the required parameters
```
nm <- as.character(sapply(colnames(df)[-(1:2)],function(v) {str_split(v,"\\.")[[1]][6]}))
nm <- as.character(sapply(nm,function(v){str_split(v,"_part")[[1]][1]}))
```
The code above is to keep the filenames of different parts consistent. This is the only part for users to change according the actual filenames.
