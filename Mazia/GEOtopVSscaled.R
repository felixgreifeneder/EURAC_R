#geotop vs stations

source("./Mazia/extr_station_SMC.R")
source("./Mazia/extr_scaled_SMC_8st.R")
source("./Mazia/extr_GEOtop_SMAPgrid.R")

plotSTATIONSvsGEOTOP <- function(){
  
  load("./Mazia/meanGEOtop2.dat")
  
  scaled_stations <- extrScaledSMC("X:/Workspaces/GrF/Processing/SMAP/scaling/8st_reprocessing_20150318")
  inSitu_stations <- extr_station_SMC("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration")
  GEOtop <- extractGEOTOP("X:/Workspaces/GrF/Processing/SMAP/GEOTOP_ts_SMAPgrid/")
  
  inSitu_NaN <- inSitu_stations
  inSitu_NaN$p1[which(inSitu_NaN$p1 == -1)] <- NaN
  inSitu_NaN$p2[which(inSitu_NaN$p2 == -1)] <- NaN
  inSitu_NaN$p3[which(inSitu_NaN$p3 == -1)] <- NaN
  inSitu_NaN$p4[which(inSitu_NaN$p4 == -1)] <- NaN
  inSitu_NaN$p5[which(inSitu_NaN$p5 == -1)] <- NaN
  inSitu_NaN$p6[which(inSitu_NaN$p6 == -1)] <- NaN
  inSitu_NaN$p7[which(inSitu_NaN$p7 == -1)] <- NaN
  inSitu_NaN$p8[which(inSitu_NaN$p8 == -1)] <- NaN
  inSitu_NaN$p9[which(inSitu_NaN$p9 == -1)] <- NaN
  inSitu_NaN$p10[which(inSitu_NaN$p10 == -1)] <- NaN
  scaled_NaN <- scaled_stations
  scaled_NaN$p1[which(scaled_NaN$p1 == -1)] <- NaN
  scaled_NaN$p2[which(scaled_NaN$p2 == -1)] <- NaN
  scaled_NaN$p3[which(scaled_NaN$p3 == -1)] <- NaN
  scaled_NaN$p4[which(scaled_NaN$p4 == -1)] <- NaN
  scaled_NaN$p5[which(scaled_NaN$p5 == -1)] <- NaN
  scaled_NaN$p6[which(scaled_NaN$p6 == -1)] <- NaN
  scaled_NaN$p7[which(scaled_NaN$p7 == -1)] <- NaN
  scaled_NaN$p8[which(scaled_NaN$p8 == -1)] <- NaN
  scaled_NaN$p9[which(scaled_NaN$p9 == -1)] <- NaN
  scaled_NaN$p10[which(scaled_NaN$p10 == -1)] <- NaN
  GEOtop_NaN <- GEOtop
  GEOtop_NaN$P1[which(GEOtop_NaN$P1 < 0.12)] <- NaN
  GEOtop_NaN$P2[which(GEOtop_NaN$P2 < 0.12)] <- NaN
  GEOtop_NaN$P3[which(GEOtop_NaN$P3 < 0.12)] <- NaN
  GEOtop_NaN$P4[which(GEOtop_NaN$P4 < 0.12)] <- NaN
  GEOtop_NaN$P5[which(GEOtop_NaN$P5 < 0.12)] <- NaN
  GEOtop_NaN$P6[which(GEOtop_NaN$P6 < 0.12)] <- NaN
  GEOtop_NaN$P7[which(GEOtop_NaN$P7 < 0.12)] <- NaN
  GEOtop_NaN$P8[which(GEOtop_NaN$P8 < 0.12)] <- NaN
  GEOtop_NaN$P9[which(GEOtop_NaN$P9 < 0.12)] <- NaN
  GEOtop_NaN$P10[which(GEOtop_NaN$P10 < 0.12)] <- NaN
  
  
  mean_insitu <- data.frame(date=inSitu_NaN$date, smc=rowMeans(inSitu_NaN[,2:11], na.rm=T))
  #mean_insitu <- data.frame(data=inSitu_NaN$date, smc=inSitu_NaN$m4)
  mean_scaled <- data.frame(date=scaled_NaN$date, smc=rowMeans(scaled_NaN[,2:11], na.rm=T))
  #mean_scaled <- data.frame(date=scaled_NaN$date, smc=scaled_NaN$p8)
  mean_GEOtop <- data.frame(date=GEOtop_NaN$date, smc=rowMeans(GEOtop_NaN[,2:11], na.rm=T))
  mean_GEOtop <- GEOtopTS
  
  
  minDate <- min(min(scaled_stations$date), min(inSitu_stations$date), min(GEOtop$date))
  maxDate <- max(max(scaled_stations$date), max(inSitu_stations$date), max(GEOtop$date))
  
  #nDays <- maxDate - minDate + 1
  DateList <- c(seq(as.Date("2010-05-01"), as.Date("2010-10-31"), by="days"), 
                seq(as.Date("2011-05-01"), as.Date("2011-10-31"), by="days"),
                seq(as.Date("2012-05-01"), as.Date("2012-10-31"), by="days"),
                seq(as.Date("2013-05-01"), as.Date("2013-10-31"), by="days"),
                seq(as.Date("2014-05-01"), as.Date("2014-10-31"), by="days"))
  nDays <- length(DateList)
  
  combinedDF <- data.frame(date=rep(as.Date(minDate), 20*nDays),
                           GEOtop=rep(-1, 20*nDays),
                           station=rep(-1, 20*nDays),
                           pixel=rep(-1, 20*nDays),
                           scaled=rep(factor("no", levels=c("no","yes")), 20*nDays))
  combinedDF_mean <- data.frame(date=rep(as.Date(minDate), 2*nDays),
                           GEOtop=rep(-1, 2*nDays),
                           station=rep(-1, 2*nDays),
                           scaled=rep(factor("no", levels=c("no","yes")), 2*nDays))
  
#--------------------------------------------------------------------------------------------------
  #STATION/PIXEL WISE
#--------------------------------------------------------------------------------------------------
  
  cntr <- 1
  
  #PSind <- c(2, 0, 1, 3, 7, 0, 0, 6, 4, 0)
  
  for (i in 1:nDays){
    
    #for (j in c(1,3,4,5,8,9)){
    for (j in c(1,2,3,4,5,6,7,8,9,10)){
      
      if (any(scaled_stations$date == DateList[i]) & any(inSitu_stations$date == DateList[i]) & any(GEOtop$date == DateList[i])){
        indInSitu <- which(inSitu_stations$date == DateList[i])
        indScaled <- which(scaled_stations$date == DateList[i])
        indGEOtop <- which(GEOtop$date == DateList[i])
        
        if (scaled_stations[indScaled, j+1] != -1 & inSitu_stations[indInSitu, j+1] != -1 & GEOtop[indGEOtop, j+1] != -9999){
          combinedDF$date[cntr] <- DateList[i]
          combinedDF$GEOtop[cntr] <- GEOtop[indGEOtop, j+1]
          combinedDF$station[cntr] <- inSitu_stations[indInSitu, j+1]
          combinedDF$pixel[cntr] <- j
          cntr <- cntr+1
          
          combinedDF$date[cntr] <- DateList[i]
          combinedDF$GEOtop[cntr] <- GEOtop[indGEOtop, j+1]
          combinedDF$station[cntr] <- scaled_stations[indScaled, j+1]
          combinedDF$pixel[cntr] <- j
          combinedDF$scaled[cntr] <- "yes"
          cntr <- cntr+1
          
        }
      }
      
    }
    
  }
  
#--------------------------------------------------------------------------------------------------
  #MEAN
#--------------------------------------------------------------------------------------------------
  
  cntr <- 1
  
  for (i in 1:nDays){
  
    if (any(mean_insitu$date == DateList[i]) & any(mean_scaled$date == DateList[i]) & any(mean_GEOtop$date == DateList[i])){
      indInSitu <- which(mean_insitu$date == DateList[i])
      indScaled <- which(mean_scaled$date == DateList[i])
      indGEOtop <- which(mean_GEOtop$date == DateList[i])
      
      if (is.finite(mean_scaled$smc[indScaled]) & is.finite(mean_insitu$smc[indInSitu]) & is.finite(mean_GEOtop$smc[indGEOtop])){
        combinedDF_mean$date[cntr] <- DateList[i]
        combinedDF_mean$GEOtop[cntr] <- mean_GEOtop$smc[indGEOtop]
        combinedDF_mean$station[cntr] <- mean_insitu$smc[indInSitu]
        cntr <- cntr+1
          
        combinedDF_mean$date[cntr] <- DateList[i]
        combinedDF_mean$GEOtop[cntr] <- mean_GEOtop$smc[indGEOtop]
        combinedDF_mean$station[cntr] <- mean_scaled$smc[indScaled]
        combinedDF_mean$scaled[cntr] <- "yes"
        cntr <- cntr+1
          
      }
    }
    
    
  }
#   
#--------------------------------------------------------------------------------------------------
  #PLOT
#--------------------------------------------------------------------------------------------------
  
  valid <- which(combinedDF$GEOtop != -1 & combinedDF$station != -1)
  print(valid)
  no_scaling <- which(combinedDF$GEOtop != -1 & combinedDF$station != -1 & combinedDF$scaled=="no")
  yes_scaling <- which(combinedDF$GEOtop != -1 & combinedDF$station != -1 & combinedDF$scaled=="yes" )
  print("no scaling:")
  print(cor(combinedDF$GEOtop[no_scaling], combinedDF$station[no_scaling]))
  print(length(no_scaling))
  print("yes scaling:")
  print(cor(combinedDF$GEOtop[yes_scaling], combinedDF$station[yes_scaling]))
  print(length(yes_scaling))
  
  ggplot(combinedDF[valid,], aes(y=GEOtop, x=station, color=scaled)) + 
    geom_point(shape="*", size=8) + 
    facet_grid(.~ pixel) +
    scale_colour_hue(l=50, name="Scaled") + 
    theme(aspect.ratio=1) +
    theme(axis.title.x = element_text(face="bold", size=18),
          axis.text.x = element_text(size=14),
          axis.title.y = element_text(face="bold", size=18),
          axis.text.y = element_text(size=14),
          strip.text.x = element_text(size=14, face="bold"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14)) +
    scale_x_continuous(name="\nSMC [m3/m3]", limits=c(0,0.4)) +
    scale_y_continuous(name="SMC [m3/m3]\n", limits=c(0,0.4)) +
    geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1)
  
  
#   valid <- which(combinedDF_mean$GEOtop != -1 & combinedDF_mean$station != -1)
#   print(valid)
#   no_scaling <- which(combinedDF_mean$GEOtop != -1 & combinedDF_mean$station != -1 & combinedDF_mean$scaled=="no")
#   yes_scaling <- which(combinedDF_mean$GEOtop != -1 & combinedDF_mean$station != -1 & combinedDF_mean$scaled=="yes" )
#   print("no scaling:")
#   print(cor(combinedDF_mean$GEOtop[no_scaling], combinedDF_mean$station[no_scaling]))
#   print(length(no_scaling))
#   print("yes scaling:")
#   print(cor(combinedDF_mean$GEOtop[yes_scaling], combinedDF_mean$station[yes_scaling]))
#   print(length(yes_scaling))
#   
#   ggplot(combinedDF_mean[valid,], aes(x=GEOtop, y=station, color=scaled)) + 
#     geom_point(shape="*", size=8) + 
#     facet_grid(.~ scaled) +
#     scale_colour_hue(l=50, name="Scaled") + 
#     theme(aspect.ratio=1) +
#     theme(axis.title.x = element_text(face="bold", size=18),
#           axis.text.x = element_text(size=14),
#           axis.title.y = element_text(face="bold", size=18),
#           axis.text.y = element_text(size=14),
#           strip.text.x = element_text(size=14, face="bold"),
#           legend.title = element_text(size=14),
#           legend.text = element_text(size=14)) +
#     scale_x_continuous(name="\nSMC [m3/m3]", limits=c(0,0.4)) +
#     scale_y_continuous(name="SMC [m3/m3]\n", limits=c(0,0.4)) 
  
  
  
}