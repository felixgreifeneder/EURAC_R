#extract a dataframe holding the station smc

extr_station_SMC <- function(path){
  
  #--------------------------------------------------------------------------------------------------
  #extract the in-situ time-series
  #--------------------------------------------------------------------------------------------------
  
  b1 <- read.table(paste(path,"/B1_2009-2014_averages_CALIBRATED.csv",sep=""), 
                   header=T, sep=",", stringsAsFactors=F)
  b2 <- read.table(paste(path,"/B2_2009-2014_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  m1 <- read.table(paste(path,"/M1_2009-2014_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  m2 <- read.table(paste(path,"/M2_2010-2013_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  m3 <- read.table(paste(path,"/M3_2009-2014_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  m4 <- read.table(paste(path,"/M4_2009-2014_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  m5 <- read.table(paste(path,"/M5_2009-2014_averages_CALIBRATED.csv",sep=""),
                   header=T, sep=",", stringsAsFactors=F)
  s4 <- read.table(paste(path,"/S4_2009-2014_averages_CALIBRATED.csv",sep=""),
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
  
  inSitu <- data.frame(date=c(inSituST:inSituET),
                       p1=rep(-1,inSituET-inSituST+1),
                       p2=rep(-1,inSituET-inSituST+1),
                       p3=rep(-1,inSituET-inSituST+1),
                       p4=rep(-1,inSituET-inSituST+1),
                       p5=rep(-1,inSituET-inSituST+1),
                       p6=rep(-1,inSituET-inSituST+1),
                       p7=rep(-1,inSituET-inSituST+1),
                       p8=rep(-1,inSituET-inSituST+1),
                       p9=rep(-1,inSituET-inSituST+1),
                       p10=rep(-1,inSituET-inSituST+1))
  
  for (dInd in inSituST:inSituET){
    
    rowIndex <- dInd-inSituST+1
    
    if (any(b2$DATE==dInd)) inSitu$p1[rowIndex] <- b2$AVGSWC5CM_CAL[which(b2$DATE==dInd)]
    if (any(s4$DATE==dInd)) inSitu$p2[rowIndex] <- s4$AVGSWC5CM_CAL[which(s4$DATE==dInd)]
    if (any(b1$DATE==dInd)) inSitu$p3[rowIndex] <- b1$AVGSWC5CM_CAL[which(b1$DATE==dInd)]
    if (any(m1$DATE==dInd)) inSitu$p4[rowIndex] <- m1$AVGSWC5CM_CAL[which(m1$DATE==dInd)]
    if (any(m5$DATE==dInd)) inSitu$p5[rowIndex] <- m5$AVGSWC5CM_CAL[which(m5$DATE==dInd)]
    if (any(m1$DATE==dInd)) inSitu$p6[rowIndex] <- m1$AVGSWC5CM_CAL[which(m1$DATE==dInd)]
    if (any(m1$DATE==dInd)) inSitu$p7[rowIndex] <- m1$AVGSWC5CM_CAL[which(m1$DATE==dInd)]
    if (any(m4$DATE==dInd)) inSitu$p8[rowIndex] <- m4$AVGSWC5CM_CAL[which(m4$DATE==dInd)]
    if (any(m2$DATE==dInd)) inSitu$p9[rowIndex] <- m2$AVGSWC5CM_CAL[which(m2$DATE==dInd)]
    if (any(m3$DATE==dInd)) inSitu$p10[rowIndex] <- m3$AVGSWC5CM_CAL[which(m3$DATE==dInd)]
    
  }
  

  inSitu$date <- as.Date(JD(inSitu$date, inverse=T))
  
  return(inSitu)
  
}