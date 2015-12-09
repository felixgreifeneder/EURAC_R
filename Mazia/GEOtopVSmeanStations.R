
library("insol")
library(ggplot2)
source("~/R/ASCAT/calculate_weights.R")


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
print(JD(inSitu$date[val[313]], inverse=T))

#create mean SMC times-series calculated over all stations
meanInSitu <- data.frame(date=as.Date(JD(inSitu$date[val], inverse=T)), SMC5CM=rowMeans(cbind(inSitu$b1[val], inSitu$b2[val], inSitu$m1[val], inSitu$m2[val],
                                                                                                                               inSitu$m3[val], inSitu$m4[val], inSitu$m5[val], inSitu$s4[val])))

#create weighted average SMC time-series
scaling_weights <- calcWeights("X:/Workspaces/GrF/01_Data/SMAP/CLASSIFICATION_8st/maxlike_nothresh.dat")
weighted_meanInSitu <- data.frame(date=as.Date(JD(inSitu$date[val], inverse=T)), SMC5CM=apply(cbind(inSitu$b1[val], inSitu$b2[val], inSitu$m1[val], inSitu$m2[val],
                                                                                                                                     inSitu$m3[val], inSitu$m4[val], inSitu$m5[val], inSitu$s4[val]),
                                                                                                                               MARGIN=1,
                                                                                                                               FUN=weighted.mean, 
                                                                                                                               w=as.vector(scaling_weights)))
rm(b1,b2,m1,m2,m3,m4,m5,s4,scaling_weights,val,rowIndex,dInd,inSituET,inSituST)

#--------------------------------------------------------------------------------------------------
#load GEOtop
#--------------------------------------------------------------------------------------------------

load("./Mazia/meanGEOtop.dat")

all <- data.frame(dat=rep(-1, nrow(meanInSitu)*2), gtop=rep(-1, nrow(meanInSitu)*2),
                  nstu=rep(-1, nrow(meanInSitu)*2), 
                  scaled=rep(factor("no", levels=c("no", "yes")), nrow(meanInSitu)*2))

cntr <- 1
for (i in 1:nrow(meanInSitu)){
  
  if (any(GEOtopTS$date == meanInSitu$date[i])){
    
    all$dat[cntr] <- meanInSitu$date[i]
    all$gtop[cntr] <- GEOtopTS$smc[which(GEOtopTS$date == meanInSitu$date[i])]
    all$nstu[cntr] <- meanInSitu$SMC5CM[i]
    all$scaled[cntr] = "no"
    cntr <- cntr + 1
    
    all$dat[cntr] <- meanInSitu$date[i]
    all$gtop[cntr] <- GEOtopTS$smc[which(GEOtopTS$date == meanInSitu$date[i])]
    all$nstu[cntr] <- weighted_meanInSitu$SMC5CM[i]
    all$scaled[cntr] = "yes"
    cntr <- cntr + 1
  }
  
}

dev.new()
ggplot(all, aes(x=gtop, y=nstu, color=scaled)) + geom_point(shape=1) + facet_grid(.~ scaled) + xlim(0,0.4) + ylim(0,0.4)
