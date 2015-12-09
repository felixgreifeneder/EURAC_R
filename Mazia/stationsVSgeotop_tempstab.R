#plot scaled and un-scaled in-situ smc data versus the mean GEOtop valuse -> temporal stability

library("insol")
library(ggplot2)


load("./Mazia/meanGEOtop2.dat")

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

rm(b1,b2,m1,m2,m3,m4,m5,s4)

val <- which(inSitu$b1 != -1 & inSitu$b2 != -1 & inSitu$m1 != -1 &
               inSitu$m2 != -1 & inSitu$m3 != -1 & inSitu$m4 != -1 &
               inSitu$m5 != -1 & inSitu$s4 != -1)

#--------------------------------------------------------------------------------------------------
#LOAD SCALED TIME-SERIES (SMAP PIXEL FOOTPRINTS)
#--------------------------------------------------------------------------------------------------

scaled_SMC <- list()

for (i in 1:10){
  
  scaled_SMC[[i]] <- read.table(paste("X:/Workspaces/GrF/Processing/SMAP/scaling/8st_reprocessing_20150318/vwc_scaled_SMAP",sprintf("%02d", i),".txt", sep=""),
                                sep=",", header=T, stringsAsFactors=F, quote="")
  
}

startDay <- min(min(scaled_SMC[[1]][,1]), 
                min(scaled_SMC[[2]][,1]),
                min(scaled_SMC[[3]][,1]),
                min(scaled_SMC[[4]][,1]),
                min(scaled_SMC[[5]][,1]),
                min(scaled_SMC[[6]][,1]),
                min(scaled_SMC[[7]][,1]),
                min(scaled_SMC[[8]][,1]),
                min(scaled_SMC[[9]][,1]),
                min(scaled_SMC[[10]][,1]))

endDay <- max(max(scaled_SMC[[1]][,1]), 
              max(scaled_SMC[[2]][,1]),
              max(scaled_SMC[[3]][,1]),
              max(scaled_SMC[[4]][,1]),
              max(scaled_SMC[[5]][,1]),
              max(scaled_SMC[[6]][,1]),
              max(scaled_SMC[[7]][,1]),
              max(scaled_SMC[[8]][,1]),
              max(scaled_SMC[[9]][,1]),
              max(scaled_SMC[[10]][,1]))

inSitu_scaled <- data.frame(date=c(startDay:endDay),
                            b1=rep(-1,endDay-startDay+1),
                            b2=rep(-1,endDay-startDay+1),
                            m1=rep(-1,endDay-startDay+1),
                            m2=rep(-1,endDay-startDay+1),
                            m3m4=rep(-1,endDay-startDay+1),
                            m5=rep(-1,endDay-startDay+1))

for (i in 1:(endDay-startDay+1)){
  
  if (any(scaled_SMC[[1]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$b2[i] <- scaled_SMC[[1]][which(scaled_SMC[[1]][,1] == inSitu_scaled$date[i]),2] 
  }
  if (any(scaled_SMC[[3]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$b1[i] <- scaled_SMC[[3]][which(scaled_SMC[[3]][,1] == inSitu_scaled$date[i]),2] 
  }
  if (any(scaled_SMC[[5]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$m5[i] <- scaled_SMC[[5]][which(scaled_SMC[[5]][,1] == inSitu_scaled$date[i]),2] 
  }
  if (any(scaled_SMC[[4]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$m1[i] <- scaled_SMC[[4]][which(scaled_SMC[[4]][,1] == inSitu_scaled$date[i]),2] 
  }
  if (any(scaled_SMC[[9]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$m2[i] <- scaled_SMC[[9]][which(scaled_SMC[[9]][,1] == inSitu_scaled$date[i]),2] 
  }
  if (any(scaled_SMC[[8]][,1] == inSitu_scaled$date[i])) { 
    inSitu_scaled$m3m4[i] <- scaled_SMC[[8]][which(scaled_SMC[[8]][,1] == inSitu_scaled$date[i]),2] 
  }
}

rm(scaled_SMC)

#--------------------------------------------------------------------------------------------------
#COMBINE AND PLOT
#--------------------------------------------------------------------------------------------------

inSitu$date <- as.Date(JD(inSitu$date, inverse=T))
inSitu_scaled$date <- as.Date(JD(inSitu_scaled$date, inverse=T))

startDay <- min(c(inSitu$date[1], inSitu_scaled$date[1], GEOtopTS$date[1]))
endDay <- max(c(inSitu$date[1717], inSitu_scaled$date[1158], GEOtopTS$date[1492]))

ext_inSitu <- data.frame(date=c(startDay:endDay),
                         b1=rep(-1,endDay-startDay+1),
                         b2=rep(-1,endDay-startDay+1),
                         m1=rep(-1,endDay-startDay+1),
                         m2=rep(-1,endDay-startDay+1),
                         m3=rep(-1,endDay-startDay+1),
                         m4=rep(-1,endDay-startDay+1),
                         m5=rep(-1,endDay-startDay+1),
                         s4=rep(-1,endDay-startDay+1))

ext_inSitu_scaled <- data.frame(date=c(startDay:endDay),
                                b1=rep(-1,endDay-startDay+1),
                                b2=rep(-1,endDay-startDay+1),
                                m1=rep(-1,endDay-startDay+1),
                                m2=rep(-1,endDay-startDay+1),
                                m3m4=rep(-1,endDay-startDay+1),
                                m5=rep(-1,endDay-startDay+1))

ext_GEOtopTS <- data.frame(date=c(startDay:endDay),
                           smc=rep(-1,endDay-startDay+1))

for (i in 1:(endDay-startDay+1)){
  if (any(GEOtopTS$date==ext_GEOtopTS$date[i])) { 
    ind <- which(GEOtopTS$date==ext_GEOtopTS$date[i]) 
    ext_GEOtopTS$smc[i] <- GEOtopTS$smc[ind]
  }
  
  if (any(inSitu$date==ext_inSitu$date[i])){
    ind <- which(inSitu$date==ext_inSitu$date[i])
    ext_inSitu[i,2:9] <- inSitu[ind,2:9]
  }
  
  if (any(inSitu_scaled$date==ext_inSitu_scaled$date[i])){
    ind <- which(inSitu_scaled$date==ext_inSitu_scaled$date[i])
    ext_inSitu_scaled[i,2:7] <- inSitu_scaled[ind,2:7]
  }
}

plotDateInd <- which(((ext_GEOtopTS$date >= as.Date("2009-05-01") & ext_GEOtopTS$date <= as.Date("2009-10-31")) |
                      (ext_GEOtopTS$date >= as.Date("2010-05-01") & ext_GEOtopTS$date <= as.Date("2010-10-31")) |
                      (ext_GEOtopTS$date >= as.Date("2011-05-01") & ext_GEOtopTS$date <= as.Date("2011-10-31")) |
                      (ext_GEOtopTS$date >= as.Date("2012-05-01") & ext_GEOtopTS$date <= as.Date("2012-10-31")) |
                      (ext_GEOtopTS$date >= as.Date("2013-05-01") & ext_GEOtopTS$date <= as.Date("2013-10-31")) |
                      (ext_GEOtopTS$date >= as.Date("2014-05-01") & ext_GEOtopTS$date <= as.Date("2014-10-31"))))

ext_GEOtopTS <- ext_GEOtopTS[plotDateInd,]
ext_inSitu <- ext_inSitu[plotDateInd,]
ext_inSitu_scaled <- ext_inSitu_scaled[plotDateInd,]

rm(plotDateInd)

all <- data.frame(date=rep(-1,811*14), 
                  model=rep(-1,811*14), 
                  insitu=rep(-1,811*14), 
                  scaled=rep(factor("no", levels = c("yes", "no")),811*14))

cntr <- 1
#length of TS'
for (i in 1:(811)){
  #unscaled
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$b1[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$b1[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$b2[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$b2[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$m1[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$m1[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$m2[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$m2[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$m3[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$m3[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$m4[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$m4[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$m5[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$m5[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu$s4[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu$s4[i]
    all$scaled[cntr] <- "no"
    cntr <- cntr + 1
  }
  
  
  #scaled
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$b1[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$b1[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$b2[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$b2[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$m1[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$m1[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$m2[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$m2[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$m3m4[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$m3m4[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
  if (ext_GEOtopTS$smc[i] != -1 & ext_inSitu_scaled$m5[i] != -1) { 
    all$date[cntr] <- ext_GEOtopTS$date[i] 
    all$model[cntr] <- ext_GEOtopTS$smc[i]
    all$insitu[cntr] <- ext_inSitu_scaled$m5[i]
    all$scaled[cntr] <- "yes"
    cntr <- cntr + 1
  }
}

dev.new()
ggplot(all, aes(x=model, y=insitu, color=scaled)) + geom_point(shape=1) + facet_grid(.~ scaled) + xlim(0,0.5) + ylim(0,0.5) + theme(aspect.ratio=1)
