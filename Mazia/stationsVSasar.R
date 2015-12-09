#asar vs stations

library(ggplot2)

source("./Mazia/extr_scaled_SMC_8st.R")
source("./Mazia/extr_ASAR_SIG0_SMAPgrid.R")
source("./Mazia/extr_station_SMC.R")

plotSTATIONSvsASAR <- function(){
  asarSIG0 <- extractASAR("X:/Workspaces/GrF/Processing/SMAP/ASAR_ts_SMAPgrid/")
  scaled_stations <- extrScaledSMC("X:/Workspaces/GrF/Processing/SMAP/scaling/8st_reprocessing_20150318")
  inSitu_stations <- extr_station_SMC("X:/Workspaces/GrF/01_Data/ANCILLARY/Station_Calibration")
  
#   asarSIG0_NaN <- asarSIG0
#   asarSIG0_NaN$SIG0_1[which(asarSIG0_NaN$SIG0_1 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_2[which(asarSIG0_NaN$SIG0_2 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_3[which(asarSIG0_NaN$SIG0_3 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_4[which(asarSIG0_NaN$SIG0_4 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_5[which(asarSIG0_NaN$SIG0_5 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_6[which(asarSIG0_NaN$SIG0_6 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_7[which(asarSIG0_NaN$SIG0_7 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_8[which(asarSIG0_NaN$SIG0_8 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_9[which(asarSIG0_NaN$SIG0_9 == -32768)] <- NaN
#   asarSIG0_NaN$SIG0_10[which(asarSIG0_NaN$SIG0_10 == -32768)] <- NaN
#   inSitu_NaN <- inSitu_stations
#   inSitu_NaN$b1[which(inSitu_NaN$b1 == -1)] <- NaN
#   inSitu_NaN$b2[which(inSitu_NaN$b2 == -1)] <- NaN
#   inSitu_NaN$m1[which(inSitu_NaN$m1 == -1)] <- NaN
#   inSitu_NaN$m2[which(inSitu_NaN$m2 == -1)] <- NaN
#   inSitu_NaN$m3[which(inSitu_NaN$m3 == -1)] <- NaN
#   inSitu_NaN$m4[which(inSitu_NaN$m4 == -1)] <- NaN
#   inSitu_NaN$m5[which(inSitu_NaN$m5 == -1)] <- NaN
#   scaled_NaN <- scaled_stations
#   scaled_NaN$p1[which(scaled_NaN$p1 == -1)] <- NaN
#   scaled_NaN$p2[which(scaled_NaN$p2 == -1)] <- NaN
#   scaled_NaN$p3[which(scaled_NaN$p3 == -1)] <- NaN
#   scaled_NaN$p4[which(scaled_NaN$p4 == -1)] <- NaN
#   scaled_NaN$p5[which(scaled_NaN$p5 == -1)] <- NaN
#   scaled_NaN$p6[which(scaled_NaN$p6 == -1)] <- NaN
#   scaled_NaN$p7[which(scaled_NaN$p7 == -1)] <- NaN
#   scaled_NaN$p8[which(scaled_NaN$p8 == -1)] <- NaN
#   scaled_NaN$p9[which(scaled_NaN$p9 == -1)] <- NaN
#   scaled_NaN$p10[which(scaled_NaN$p10 == -1)] <- NaN
#   
#   mean_insitu <- data.frame(date=inSitu_NaN$date, smc=rowMeans(inSitu_NaN[,2:9], na.rm=T))
#   mean_scaled <- data.frame(date=scaled_NaN$date, smc=rowMeans(scaled_NaN[,2:9], na.rm=T))
#   asarSIG0_NaN[,2:11] <- 10^(asarSIG0_NaN[,2:11]/10)
#   mean_asar <- data.frame(date=asarSIG0_NaN$Day, sig0=rowMeans(asarSIG0_NaN[,2:11], na.rm=T))
#   mean_asar[,2] <- 10*log10(mean_asar[,2])
  
  minDate <- min(min(scaled_stations$date), min(inSitu_stations$date), min(asarSIG0$Day))
  maxDate <- max(max(scaled_stations$date), max(inSitu_stations$date), max(asarSIG0$Day))
  minDate <- as.Date("2010-05-01")
  in1 <- as.Date("2010-10-31")
  in2 <- as.Date("2011-05-01")
  maxDate <- as.Date("2011-10-31")
  
  #nDays <- maxDate - minDate + 1
  DateList <- c(seq(minDate, in1, by="days"), seq(in2, maxDate, by="days"))
  nDays <- length(DateList)
  
  
  sumry <- data.frame(date=rep(as.Date(minDate),18*nDays),
                        smc=rep(-1, 18*nDays),
                        sig0=rep(-1, 18*nDays),
                        pixel=rep(0, 18*nDays),
                        scaled=rep(factor("no", levels=c("no","yes")),18*nDays))
  
  sumry_mean <- data.frame(date=rep(as.Date(minDate),2*nDays),
                           sig0=rep(-1, 2*nDays),
                           smc=rep(-1, 2*nDays),
                           scaled=rep(factor("yes", levels=c("no","yes")), 2*nDays))

  #------------------------------------------------------------------------------------------------
  #PIXEL/STATION - WISE
  #------------------------------------------------------------------------------------------------
    
  cntr <- 1
  
  #iDinList <- c(2,0,1,3,7,0,0,5,4,0)
  
  for (i in 1:nDays){
    
    for (j in c(1,2,3,7,8,9,10)){
      
      if (any(asarSIG0$Day == DateList[i]) & any(scaled_stations$date == DateList[i]) & any(inSitu_stations$date == DateList[i])){
        
        idA <- which(asarSIG0$Day == DateList[i])
        idS <- which(scaled_stations$date == DateList[i])
        idI <- which(inSitu_stations$date == DateList[i])
        
        if (asarSIG0[idA, j+1] != -32768 & scaled_stations[idS, j+1] !=-1 & inSitu_stations[idI, j+1] != -1){
          sumry$date[cntr] <- DateList[i]
          sumry$sig0[cntr] <- asarSIG0[idA, j+1]
          sumry$smc[cntr] <- inSitu_stations[idI, j+1]
          sumry$pixel[cntr] <- j
          cntr <- cntr + 1
          sumry$date[cntr] <- DateList[i]
          sumry$sig0[cntr] <- asarSIG0[idA, j+1]
          sumry$smc[cntr] <- scaled_stations[idS, j+1]
          sumry$pixel[cntr] <- j
          sumry$scaled[cntr] <- "yes"
          cntr <- cntr + 1
        }
        
#         if (asarSIG0[idA, j+1] != -32768 & inSitu_stations[idI, j+1] != -1) {
#           sumry$date[cntr] <- DateList[i]
#           sumry$sig0[cntr] <- asarSIG0[idA, j+1]
#           sumry$smc[cntr] <- inSitu_stations[idI, j+1]
#           sumry$pixel[cntr] <- j
#           cntr <- cntr + 1   
#         }
        
      }
        
    }
  }

  
  #------------------------------------------------------------------------------------------------
  #MEAN
  #------------------------------------------------------------------------------------------------
  
#   for (i in 1:nDays){
#     if (any(mean_insitu$date == sumry_mean$date[i]) & any(mean_scaled$date == sumry_mean$date[i]) & any(mean_asar$date == sumry_mean$date[i])){
#       
#       asar_ind <- which(mean_asar$date == sumry_mean$date[i])
#       insitu_ind <- which(mean_insitu$date == sumry_mean$date[i])
#       scaled_ind <- which(mean_scaled$date == sumry_mean$date[i])
#       
#       if (is.finite(mean_asar$sig0[asar_ind]) & is.finite(mean_insitu$smc[insitu_ind]) & is.finite(mean_scaled$smc[scaled_ind])){
#         
#         sumry_mean$sig0[i] <- mean_asar$sig0[asar_ind]
#         sumry_mean$smc[i] <- mean_insitu$smc[insitu_ind]
#         sumry_mean$scaled[i] <- "no"
#         sumry_mean$sig0[nDays+i] <- mean_asar$sig0[asar_ind]
#         sumry_mean$smc[nDays+i] <- mean_scaled$smc[scaled_ind]
#         
#       }
#       
#     }
#   }
#   
#--------------------------------------------------------------------------------------------------
  #PLOTTING
#--------------------------------------------------------------------------------------------------
  
  valid <- which(sumry$sig0 != -1 & sumry$smc != -1 & sumry$sig0 < -4)
  no_scaling <- which(sumry$sig0 != -1 & sumry$smc != -1 & sumry$scaled=="no" & sumry$sig0 < -4)
  yes_scaling <- which(sumry$sig0 != -1 & sumry$smc != -1 & sumry$scaled=="yes" & sumry$sig0 < -4)
  print("no scaling:")
  print(cor(sumry$sig0[no_scaling], sumry$smc[no_scaling]))
  print(length(no_scaling))
  print("yes scaling:")
  print(cor(sumry$sig0[yes_scaling], sumry$smc[yes_scaling]))
  print(length(yes_scaling))
  
#   dev.new()
  sumry$pixel <- as.factor(sumry$pixel)
  ggplot(sumry[valid,], aes(x=sig0, y=smc, color=scaled)) + 
    geom_point(shape=1, size=3) + 
    facet_grid(.~ scaled) +
    scale_colour_hue(l=50, guide=F) + 
    theme(aspect.ratio=1) +
    theme(axis.title.x = element_text(face="bold", size=18),
          axis.text.x = element_text(size=14),
          axis.title.y = element_text(face="bold", size=18),
          axis.text.y = element_text(size=14),
          strip.text.x = element_text(size=14, face="bold"),
          legend.title = element_text(size=14),
          legend.text = element_text(size=14)) +
    geom_smooth(method=glm, se=T, fullrange=T, linetype="dashed", size=1) +
    scale_x_continuous(name="\nSIG0 [dB]") +
    scale_y_continuous(name="SMC [m3/m3]\n") 
  
#   dev.new()
#   ggplot(sumry_mean, aes(x=sig0, y=smc, color=scaled)) + geom_point(shape=1) + scale_colour_hue(l=50) + 
#     theme(aspect.ratio=1) + 
#     facet_grid(.~ scaled) +
#     geom_smooth(method=glm, se=F, fullrange=T)
  
  
  
}

