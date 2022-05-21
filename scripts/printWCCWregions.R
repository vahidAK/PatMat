#Adam's WC_regions function for StrandPhaseR

library(GenomicRanges)


printWCregions <- function(datapath, file=NULL, regionSize=5000000) {

	files <- list.files(datapath, pattern=".RData$", full=T)

	ranges <- GRangesList()
	for (filename in files) {
		counts <- get(load(filename))$counts
				
		counts.filt <- counts[width(counts) >= regionSize & counts$states == 'wc']
		ranges[[filename]] <- counts.filt
		
	}	

	ranges <- unlist(ranges)
	ranges <- sort(ranges)
	filenames <- gsub(".RData", "", basename(names(ranges)))
	names(ranges) <- NULL
	ranges$filename <- filenames

	df <- as.data.frame(ranges)
	df2print <- df[,c('seqnames','start','end','filename')]
	write.table(df2print, file, sep=":", col.names=FALSE, row.names=FALSE, quote=FALSE)
}



printWCregions("./BPR_output/data", file="./WC_CW/wc_regions.txt", regionSize=100000)
