#Runs StrandPhaseR

library(StrandPhaseR)
library(BSgenome.Hsapiens.UCSC.hg38)

args <- commandArgs(T)

strandPhaseR(inputfolder = paste0(args[2],'/../'), outputfolder = paste0(args[2],"/SPR_output"), numCPU = as.integer(args[1]),
	chromosomes = 
c('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX'),
	WCregions = './wc_regions.txt',
	pairedEndReads=as.logical(args[4]),
	positions = args[3],
	num.iterations = 4,
	exportVCF = 'wc_cw',
	bsGenome = 'BSgenome.Hsapiens.UCSC.hg38')
  
