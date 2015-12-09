#derive ERA modelled SMC daily means

library(ncdf4)

meanERA_SMC <- function(path, latlon){
  
  era_file <- nc_open(path)
  
  #read data from file
  lat <- ncvar_get(era_file, varid="latitude")
  lon <- ncvar_get(era_file, varid="longitude")
  
  timeorig <- as.POSIXct("1900-01-01 00:00:00")
  time <- as.POSIXct(ncvar_get(era_file, varid="time")*60*60, tz="UTC", origin=timeorig)
  
  smc <- ncvar_get(era_file, varid="swvl1")
  
  #determin the pixel closest to stations
  dist <- matrix(nrow=nrow(smc), ncol=ncol(smc))
  for (yi in 1:length(lon)){
    for (xi in 1:length(lat)){
      dist[yi,xi] <- sqrt((lat[xi]-latlon[1])^2+(lon[yi]-latlon[2])^2)
    }
  }
  
  #extract nearest time series
  mindist <- arrayInd(which.min(dist), .dim = c(nrow(smc), ncol(smc)))
  smcts <- smc[mindist[1],mindist[2],]
  smcts[smcts == -32767] <- NA
  
  #create daily mean smc
  days <- seq(trunc(min(time), units="days"), trunc(max(time), units="days"), "days")
  smcDaily <- rep(-1, length(days))
  
  for (ti in 1:length(days)){
    
    dT <- which(time >= days[ti] & time < (days[ti]+(24*60*60)))
    smcDaily[ti] <- mean(smcts[dT], na.rm=T)
    
  }

  
  nc_close(era_file)
  
  out <- data.frame(date=as.Date(days), smc=smcDaily, stringsAsFactors = F)
  return(out)
  
}