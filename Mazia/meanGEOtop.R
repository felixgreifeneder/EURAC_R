#computes mean geotop
library("raster")

meanGEOtop <- function(fpath){
  
  GEOtop_files <- list.files(path=fpath, pattern="thetaliqL0001", full.names=T)
  
  startday <- as.Date("2010-05-01")
  endday <- startday+length(GEOtop_files)-1
  
  

  mean_smc <- data.frame(date=seq.Date(startday, endday, by="day"), smc=rep(-1, length(GEOtop_files)))
  
  for (i in 1:length(GEOtop_files)){
    
    print(i)
    GEOtop_raster <- raster(GEOtop_files[i])
    NAvalue(GEOtop_raster) <- -9999
    
    GEOtop_raster[GEOtop_raster < 0.1] <- NaN
    mean_smc$smc[i] <- cellStats(GEOtop_raster, stat="mean", na.rm=T)
    
  }
  
  return(mean_smc)
  
}