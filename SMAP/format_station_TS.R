#create station time series

SMCts <- function(filepath, stNr, outpath){
  
  #ssearch for all station files
  filelist <- list.files(path=filepath, pattern=stNr, full.names=T)
  
  combinedTable <- data.frame()
  
  for (fpath in filelist){
    
    fileContent <- read.table(fpath, sep=",", quote="", stringsAsFactors = F, skip=3)
    combinedTable <- rbind(combinedTable, fileContent)
    
  }
  
  colnames(combinedTable) <- c("ID","Yr","Mo","Dy", "Hr","Min","TOY","Tair","PREC","SM05_1","SM05_2","SM20_1","SM20_2",
                               "ST","ST","ST","ST","ST")
  
  save(combinedTable, file=paste(outpath,"FullTS",stNr,"_scaled.dat",sep=""))
  
  return("Success")
  
}