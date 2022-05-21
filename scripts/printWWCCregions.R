#modified from Adam's WC_regions function for StrandPhaseR

library(GenomicRanges)


printWWCCregions <- function(datapath, file=NULL, regionSize=5000000, state='cc') {

	files <- list.files(datapath, pattern=".RData$", full=T)

	ranges <- GRangesList()
	for (filename in files) {
		counts <- get(load(filename))$counts
				
		counts.filt <- counts[width(counts) >= regionSize & counts$states == state]
		ranges[[filename]] <- counts.filt
		
	}	

	ranges <- unlist(ranges)
	ranges <- sort(ranges)
	filenames <- gsub(".RData", "", basename(names(ranges)))
	names(ranges) <- NULL
	ranges$filename <- filenames

	df <- as.data.frame(ranges)
	df2print <- df[,c('seqnames','start','end','filename')]
	write.table(df2print, file, sep="	", col.names=FALSE, row.names=FALSE, quote=FALSE)
}



printWWCCregions("./BPR_output/data", file="./cc_regions.txt", regionSize=100000, state='cc')
printWWCCregions("./BPR_output/data", file="./ww_regions.txt", regionSize=100000, state='ww')
