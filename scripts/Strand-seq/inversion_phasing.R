library(breakpointR)

args <- commandArgs(T)
chr <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7","chr8","chr9", "chr10","chr11","chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")


pdf()
breakpointr(inputfolder="./", outputfolder="./BPR_output", pairedEndReads=as.logical(args[3]), numCPU=as.integer(args[4]),windowsize=8000000,binMethod="size", chromosomes=chr,
	background=0.1, maskRegions=paste0(args[7],"blacklist.highdepth.centromeres.bed"))
dev.off()


exportRegions("./BPR_output/data/", file="wc_regions.txt", collapseInversions=FALSE, minRegionSize=5000000, state="wc")

library(StrandPhaseR)
library(BSgenome.Hsapiens.UCSC.hg38)

strandPhaseR(inputfolder="./", outputfolder="./SPR_output/", numCPU=as.integer(args[4]), positions=args[2], WCregions="wc_regions.txt", chromosomes=chr, num.iterations=3, exportVCF=args[1], 
	bsGenome='BSgenome.Hsapiens.UCSC.hg38', splitPhasedReads=TRUE, assume.biallelic=TRUE, pairedEndReads=as.logical(args[3]))

correctInvertedRegionPhasing(outputfolder="./", inv.bed=paste0(args[5],"/",args[1],".inv.bed"), strandphaseR.data=paste0(args[5], "/SPR_output/data/"), 
	breakpointR.data=paste0(args[5],"/BPR_output/data/"), vcfs.files=paste0(args[5],"/SPR_output/VCFfiles/"), 
	snv.positions=paste0(args[2]), input.bams=args[5], chromosomes=chr, bsGenome='BSgenome.Hsapiens.UCSC.hg38', 
	ref.fasta=args[6], recall.phased = TRUE, het.genotype = 'lenient', pairedEndReads = as.logical(args[3]), min.mapq = 10, background = 0.1, lookup.bp = 1000000, assume.biallelic = TRUE, 
	lookup.blacklist=paste0(args[7],"blacklist.highdepth.centromeres.bed"))


