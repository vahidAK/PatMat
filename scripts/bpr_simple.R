#Runs BreakpointR

library(breakpointR)

args <- commandArgs(T)

chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8","chr9", "chr10","chr11","chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
chr <- chr[1:22] 


#Run breakpointr and save to PDF
#Note that the blacklist is specified here, and may need to be changed in the rare case that a very large inversion (e.g. NA19239 on chrY) interferes with file creation. 
pdf()
breakpointr(inputfolder="./", outputfolder=args[6], pairedEndReads=as.logical(args[4]), numCPU=as.integer(args[3]), 
chromosomes=chr,windowsize=as.integer(args[1]),binMethod="size", minReads=as.integer(args[2]), background=0.2,maskRegions=paste0(args[5],"blacklist.highdepth.centromeres.bed"))
dev.off()

