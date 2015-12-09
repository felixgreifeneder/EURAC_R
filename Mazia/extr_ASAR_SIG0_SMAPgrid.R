#extract time-series of ASAR SIG0, resampled to SMAP grid

library("insol")

extractASAR <- function(path){
  
  sig0_tab <- read.table(paste(path, "ASAR_SIG0_SMAP_valmask.txt", sep=""),
                         header=T,
                         sep=",",
                         quote="",
                         stringsAsFactors=F)
  
  sig0_tab$Day <- as.Date(JD(sig0_tab$Day, inverse=T))
  
  return(sig0_tab)
  
}