#this routine plots ascat smc vs station smc to investigate the potential of the scaling function

library("ncdf4")
library("insol")
source("~/R/ASCAT/calculate_weights.R")

ascatTS <- nc_open("C:/Users/FGreifeneder/Documents/tmp_proc/ASCAT/TUW_METOP_ASCAT_WARP55R22_1395.nc")

lat <- ncvar_get(ascatTS, varid="lat")
lon <- ncvar_get(ascatTS, varid="lon")
id <- ncvar_get(ascatTS, varid="location_id")
row_size <- ncvar_get(ascatTS, varid="row_size")

sm <- ncvar_get(ascatTS, varid="sm")
time <- ncvar_get(ascatTS, varid="time")

#find location of mazia grid-point
mazia_PID <- 2369207
loc <- which(id == mazia_PID)
mazia_lat <- lat[loc]
mazia_lon <- lon[loc]
                 
print(paste("Extracting data for:", mazia_PID))
print(paste("Lat:", mazia_lat))
print(paste("Lon:", mazia_lon))

#determine start of according time-seris
start <- sum(row_size[1:(loc-1)])+1
end <- sum(row_size[1:loc])

#extract data
mazia_sm_ts <- sm[start:end]
mazia_time <- time[start:end]

#extract data for given time-frame
timOrigin <- as.POSIXct("1900-01-01 00:00:00", format="%F %H:%M:%S", tz="GMT")
mazia_time <- as.POSIXct(mazia_time*24*60*60, origin=timOrigin, tz="GMT")

startTime <- as.POSIXct("2009-05-01", tz="GMT")
endTime <- as.POSIXct("2014-10-31", tz="GMT")
timeID <- which((mazia_time >= as.POSIXct("2009-05-01", tz="GMT") & mazia_time <= as.POSIXct("2009-10-31", tz="GMT")) | 
                (mazia_time >= as.POSIXct("2010-05-01", tz="GMT") & mazia_time <= as.POSIXct("2010-10-31", tz="GMT")) |
                (mazia_time >= as.POSIXct("2011-05-01", tz="GMT") & mazia_time <= as.POSIXct("2011-10-31", tz="GMT")) |
                (mazia_time >= as.POSIXct("2012-05-01", tz="GMT") & mazia_time <= as.POSIXct("2012-10-31", tz="GMT")) |
                (mazia_time >= as.POSIXct("2013-05-01", tz="GMT") & mazia_time <= as.POSIXct("2013-10-31", tz="GMT")) |
                (mazia_time >= as.POSIXct("2014-05-01", tz="GMT") & mazia_time <= as.POSIXct("2014-10-31", tz="GMT")))

#extract time series for 2012
mazia_sm_ts <- mazia_sm_ts[timeID]
mazia_time <- mazia_time[timeID]

#--------------------------------------------------------------------------------------------------
#extract the in-situ time-series
#--------------------------------------------------------------------------------------------------

b1 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/B1_2009-2014_averages_CALIBRATEDgrl.csv", 
                 header=T, sep=",", stringsAsFactors=F)
b2 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/B2_2009-2014_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
m1 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/M1_2009-2014_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
m2 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/M2_2010-2013_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
m3 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/M3_2009-2014_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
m4 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/M4_2009-2014_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
m5 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/M5_2009-2014_averages_CALIBRATEDgrl.csv",
                 header=T, sep=",", stringsAsFactors=F)
s4 <- read.table("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration/S4_2009-2014_averages_CALIBRATED.csv",
                 header=T, sep=",", stringsAsFactors=F)

inSituST <- min(c(min(b1$DATE),
                  min(b2$DATE),
                  min(m1$DATE),
                  min(m2$DATE),
                  min(m3$DATE),
                  min(m4$DATE),
                  min(m5$DATE),
                  min(s4$DATE)))
inSituET <- max(c(max(b1$DATE),
                  max(b2$DATE),
                  max(m1$DATE),
                  max(m2$DATE),
                  max(m3$DATE),
                  max(m4$DATE),
                  max(m5$DATE),
                  max(s4$DATE)))

#merge in-situ vectors
inSitu <- data.frame(date=c(inSituST:inSituET),
                     b1=rep(-1,inSituET-inSituST+1),
                     b2=rep(-1,inSituET-inSituST+1),
                     m1=rep(-1,inSituET-inSituST+1),
                     m2=rep(-1,inSituET-inSituST+1),
                     m3=rep(-1,inSituET-inSituST+1),
                     m4=rep(-1,inSituET-inSituST+1),
                     m5=rep(-1,inSituET-inSituST+1),
                     s4=rep(-1,inSituET-inSituST+1))

for (dInd in inSituST:inSituET){
  
  rowIndex <- dInd-inSituST+1
  
  if (any(b1$DATE==dInd)) inSitu$b1[rowIndex] <- b1$AVGSWC5CM_CAL[which(b1$DATE==dInd)]
  if (any(b2$DATE==dInd)) inSitu$b2[rowIndex] <- b2$AVGSWC5CM_CAL[which(b2$DATE==dInd)]
  if (any(m1$DATE==dInd)) inSitu$m1[rowIndex] <- m1$AVGSWC5CM_CAL[which(m1$DATE==dInd)]
  if (any(m2$DATE==dInd)) inSitu$m2[rowIndex] <- m2$AVGSWC5CM_CAL[which(m2$DATE==dInd)]
  if (any(m3$DATE==dInd)) inSitu$m3[rowIndex] <- m3$AVGSWC5CM_CAL[which(m3$DATE==dInd)]
  if (any(m4$DATE==dInd)) inSitu$m4[rowIndex] <- m4$AVGSWC5CM_CAL[which(m4$DATE==dInd)]
  if (any(m5$DATE==dInd)) inSitu$m5[rowIndex] <- m5$AVGSWC5CM_CAL[which(m5$DATE==dInd)]
  if (any(s4$DATE==dInd)) inSitu$s4[rowIndex] <- s4$AVGSWC5CM_CAL[which(s4$DATE==dInd)]
  
}

val <- which(inSitu$b1 != -1 & inSitu$b2 != -1 & inSitu$m1 != -1 &
             inSitu$m2 != -1 & inSitu$m3 != -1 & inSitu$m4 != -1 &
             inSitu$m5 != -1 & inSitu$s4 != -1)

print(JD(inSitu$date[val[1]], inverse=T))
print(JD(inSitu$date[length(val)], inverse=T))

#create mean SMC times-series calculated over all stations
meanInSitu <- data.frame(date=as.POSIXct(round(JD(inSitu$date[val], inverse=T), units="days"),tz="GMT"), SMC5CM=rowMeans(cbind(inSitu$b1[val], inSitu$b2[val], inSitu$m1[val], inSitu$m2[val],
                                                                          inSitu$m3[val], inSitu$m4[val], inSitu$m5[val], inSitu$s4[val])))

#create weighted average SMC time-series
scaling_weights <- calcWeights("X:/Workspaces/GrF/01_Data/SMAP/CLASSIFICATION_8st/maxlike_nothresh.dat")
weighted_meanInSitu <- data.frame(date=as.POSIXct(round(JD(inSitu$date[val], inverse=T), units="days"),tz="GMT"), SMC5CM=apply(cbind(inSitu$b1[val], inSitu$b2[val], inSitu$m1[val], inSitu$m2[val],
                                                                                           inSitu$m3[val], inSitu$m4[val], inSitu$m5[val], inSitu$s4[val]),
                                                                                     MARGIN=1,
                                                                                     FUN=weighted.mean, 
                                                                                     w=as.vector(scaling_weights)))

timeID <- which(meanInSitu$date >= startTime & meanInSitu$date <= endTime)
meanInSitu <- meanInSitu[timeID,]
weighted_meanInSitu <- weighted_meanInSitu[timeID,]




#--------------------------------------------------------------------------------------------------
#create combined table
#--------------------------------------------------------------------------------------------------

load("./Mazia/meanGEOtop.dat")

INSITU_ASCAT_GEOTOP <- data.frame(date=round(mazia_time, units="days"), insitusmc05=rep(-1, length(mazia_time)), insitusmc05_W=rep(-1, length(mazia_time)),
                                  ascatsmc=mazia_sm_ts, geotop=rep(-1, length(mazia_time)))

cntr <- 0
for (i in 1:length(mazia_time)){
  ind <- round(mazia_time[i], unit="days")
  if (any(meanInSitu$date == ind)){
    INSITU_ASCAT_GEOTOP$insitusmc05[i] <- meanInSitu$SMC5CM[which(meanInSitu$date==ind)]
    INSITU_ASCAT_GEOTOP$insitusmc05_W[i] <- weighted_meanInSitu$SMC5CM[which(weighted_meanInSitu$date==ind)]
  }
  if (any(GEOtopTS$date == as.Date(ind))){
    INSITU_ASCAT_GEOTOP$geotop[i] <- GEOtopTS$smc[which(GEOtopTS$date==as.Date(ind))]
  }
  
}

tmpval <- which(INSITU_ASCAT_GEOTOP$insitusmc05 != -1)



#--------------------------------------------------------------------------------------------------
#comparison with GEOtop
#--------------------------------------------------------------------------------------------------

ggplot(INSITU_ASCAT_GEOTOP[tmpval,], aes(x=insitusmc05, y=ascatsmc)) + geom_point(shape=1)



