library(invertyper)

args <- commandArgs(T)

invertyper(paste0(args[1],".WW_CC.bam"), paste0(args[1],".WC_CW.bam"),
bed=paste0(args[4],"sup_table_24_inversions.including_half_intervals.bed"),
blacklist=paste0(args[4],"blacklist.highdepth.centromeres.bed"), paired_reads=as.logical(args[3]), sex=args[2], confidence=0.95,prior=c(0.9866,0.0067,0.0067),
prior_male=c(0.9866, 0.0134), output_file=paste0(args[1],".catalogue.genotyped.txt"), adjust_method="all")
