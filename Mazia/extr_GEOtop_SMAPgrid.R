#extract GEOtop time-series, which were resampled to the SMAP grid

library("insol")

extractGEOTOP <- function(path){
  
  SMC_tab <- read.table(paste(path, "GEOTOP_SMC_SMAP.txt", sep=""),
                         header=F,
                         sep=",",
                         quote="",
                         stringsAsFactors=F)
  
  SMC_tab$V1 <- as.Date(JD(SMC_tab$V1, inverse=T))
  
  colnames(SMC_tab) <- c("date", "P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10")
  
  return(SMC_tab)
  
}