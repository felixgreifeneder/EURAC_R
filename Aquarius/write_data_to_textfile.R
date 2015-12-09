#write a dataframe containing either the entire dataset or the trainingset to a text file

writetxtData <- function(dataFile, outPath){
  
  write.table(dataFile, file=outPath, quote=FALSE, sep=",", row.names=F)
  
}