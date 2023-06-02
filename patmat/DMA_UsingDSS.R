#! /usr/bin/env Rscript
# coding=utf-8
# This script performs differential methylation analysis (DMA) using the DSS R package.

suppressPackageStartupMessages(library("sys"))
suppressPackageStartupMessages(library("R.utils"))
suppressPackageStartupMessages(library("DSS"))
suppressPackageStartupMessages(library("tibble"))
args <- commandArgs(trailingOnly = TRUE)
options(scipen = 15)
print(paste("Performing differential methylation analysis",
            "using DSS version:", packageVersion("DSS")))

file1= read.table(args[1],header=TRUE,sep = '\t')
file1= file1[,c(1,2,4,5)]
colnames(file1)= c("chr","pos","N","X")

file2= read.table(args[2],header=TRUE,sep = '\t')
file2= file2[,c(1,2,4,5)]
colnames(file2)= c("chr","pos","N","X")

DMLout_results= paste(args[3],"_callDML.tsv",sep = "")
test_out= paste(args[3],"_DMLtest.tsv",sep = "")

ed= args[4]
sf=args[5]
ss=as.integer(args[6])
if (ed == "FALSE"){ ed= FALSE }
if (ed == "TRUE"){ ed= TRUE }

DSObject<- makeBSseqData(list(file1,file2),c("C1", "N1"))

if (sf == "FALSE"){
  test<- DMLtest(DSObject, 
                 group1=c("C1"), 
                 group2=c("N1"), 
                 equal.disp = ed,
                 smoothing = FALSE,
                 ncores=as.integer(args[9]))
}

if (sf == "TRUE"){
  test<- DMLtest(DSObject, 
                 group1=c("C1"), 
                 group2=c("N1"), 
                 equal.disp = ed,
                 smoothing = TRUE, 
                 smoothing.span = ss,
                 ncores=as.integer(args[9]))
}

DM_loci<- callDML(test, delta= as.double(args[7]), p.threshold= as.double(args[8]))

test<- test[order(test[,1],test[,2]),]
test<- add_column(test, pos_end = test[,2]+1, .after = 2)

DM_loci<- DM_loci[order(DM_loci[,1],DM_loci[,2]),]
DM_loci<- add_column(DM_loci, pos_end = DM_loci[,2]+1, .after = 2)

write.table(test, test_out, sep="\t", row.names=F, quote=F)
write.table(DM_loci, DMLout_results, sep="\t", row.names=F, quote=F)

