#Runs BreakpointR

library(breakpointR)

args <- commandArgs(T)

#A few unecessarly elaborate lines so that we get composite files for the correct chromosomes
chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8","chr9", "chr10","chr11","chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")
if (args[5] == "male" && args[6] == "WC") {
	chr <- chr[1:22] 
} else if ( args[5] == "female" ) {
	chr <- chr[1:23]
}


if(args[8]){
	blacklist_type=".with_chr8_inv"
}else{
	blacklist_type=""
}

#Run breakpointr and save to PDF Note that the blacklist is specified here, and may need to be changed in the rare case that a very large inversion (e.g. NA19239 on chrY) interferes with file creation.
pdf()
breakpointr(inputfolder="./", outputfolder="./BPR_output", pairedEndReads=as.logical(args[4]), numCPU=as.integer(args[3]), 
chromosomes=chr,windowsize=as.integer(args[1]),binMethod="size", minReads=as.integer(args[2]), background=0.2,maskRegions=paste0(args[7],"blacklist.highdepth.centromeres",blacklist_type,".bed"))
dev.off()

