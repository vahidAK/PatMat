#Extracts the "counts" data (i.e. strand state for segments between breakpoints) from BreakpointR output

library(breakpointR)

args <- commandArgs(T)

lww <- list.files(pattern=paste0(args[1],".+WW_CC.+RData$"))
lwc <- list.files(pattern=paste0(args[1],".+WC_CW.+RData$"))
#print(lww)
wc <- loadFromFiles(lwc)
ww <- loadFromFiles(lww)

#print(head(lapply(wc,function(x)x$counts[x$counts$states!="wc"])))
write.table(lapply(wc,function(x)x$counts[x$counts$states!="wc"]), file=paste0("../../",args[1],"_",args[2],".wc.bpr.txt"),sep="\t",quote=FALSE,row.names=FALSE)

write.table(lapply(ww,function(x)x$counts[x$counts$states!="ww"]), file=paste0("../../",args[1],"_",args[2],".ww.bpr.txt"),sep="\t",quote=FALSE,row.names=FALSE)






