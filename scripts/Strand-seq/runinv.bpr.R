library(invertyper)

args <- commandArgs(T)

invertyper(paste0(args[1],".WW_CC.bam"), paste0(args[1],".WC_CW.bam"),
bed=paste0(args[1],".bpr.txt"),
blacklist=paste0(args[4],"blacklist.highdepth.centromeres.bed"), paired_reads=as.logical(args[3]), sex=args[2], confidence=0.95,prior=c(0.9,0.05,0.05),
prior_male=c(0.9,0.1), output_file=paste0(args[1],".bpr.genotyped.txt"), adjust_method="merge")
