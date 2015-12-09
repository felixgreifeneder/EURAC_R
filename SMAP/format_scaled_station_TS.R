#create station time series

SMCts_scaled <- function(filepath, stNr, outpath){
  
  #ssearch for all station files
  filelist <- list.files(path=filepath, pattern=stNr, full.names=T)
  
  combinedTable <- data.frame()
  
  for (fpath in filelist){
    
    fileContent <- read.table(fpath, sep=",", quote="", stringsAsFactors = F, skip=4)
    combinedTable <- rbind(combinedTable, fileContent)
    
  }
  
  colnames(combinedTable) <- c("ID","Yr","Mo","Dy", "Hr","Min","TOY", "WASM_1","SWASM_1","WASM_2","SWASM_2",
                               "SMB1_1","SMB2_1","SMB3_1","SMB1_2","SMB2_2","SMB3_2")
  
  save(combinedTable, file=paste(outpath,"FullTS",stNr,"_scaled.dat",sep=""))
  
  return("Success")
  
}