#extract values for validation
library("e1071")
library("ggplot2")
library("R.matlab")

#load required data
#--------------------------------------------------------------------------------------------------
val_points <- read.table("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/Datapoints_validation_RFI.csv", 
                         header=TRUE, sep=",")
# load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/testset_noRFIfilterR")
# load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/SVR_ascatM_aqu/SVR_R")
# testset <- testset_noRFIfilter
# rm(testset_noRFIfilter)

Cdata <- readMat("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/SCA_Eumetsat_Claudia/new test data noRFI/estim_smc_vv_f_vhvv_ts_n2.mat")
testset <- as.data.frame(Cdata$datats.n[,2:17])
colnames(testset) <- c("lat",
                       "lon",
                       "time",
                       "ascatF",
                       "ascatincF",
                       "ascatM",
                       "ascatincM",
                       "ascatA",
                       "ascatincA",
                       "slope",
                       "aquVV",
                       "aquVH",
                       "aquINC",
                       "aquSFLAG",
                       "smc",
                       "TEMP")

SMCpredicted <- t(Cdata$estim.smc.vv.f.vhvv.ts.n2)

#plot
# map <- borders("world", colour="gray50", fill="gray50")
# dev.new(width=10, height=6, noRStudioGD=T)
# gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=n), data=prediction_performance, size=1.2)
# gridPlt <- gridPlt + geom_point(aes(x=lon, y=lat), data=val_points, colour="red", size=2)
# gridPlt + ggtitle("Global grid\n")


#predict
#--------------------------------------------------------------------------------------------------
# tunedModel <- SVRtuning$best.model
# 
# SMCpredicted <- predict(tunedModel, testset)


#find the matching time series
#--------------------------------------------------------------------------------------------------
# locations_filename <- "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/timeseries/locations.txt"
# locations_frame <- data.frame(name=character(nrow(val_points)), st_lat=numeric(nrow(val_points)), st_lon=numeric(nrow(val_points)), aqu_lat=numeric(nrow(val_points)), aqu_lon=numeric(nrow(val_points)))

for (station_ind in c(1:nrow(val_points))){
  
  dist_list <- numeric(0)
  
  uniqueLocations <- unique(testset[,1:2])
  nUL <- nrow(uniqueLocations)
  
  #find the closest station
  for (grid_ind in c(1:nUL)){
    
    gridLat <- uniqueLocations$lat[grid_ind]
    gridLon <- uniqueLocations$lon[grid_ind]
    stationLat <- val_points$lat[station_ind]
    stationLon <- val_points$lon[station_ind]
    
    dist <- sqrt((stationLat-gridLat)^2 + (stationLon-gridLon)^2)
    dist_list <- c(dist_list, dist)
  }
  
  mindist.ind <- which.min(dist_list)
  mindist <- min(dist_list)
  print(paste(val_points$station[station_ind], mindist))
  
  if (mindist <= 1){
    smcts.ind <- which(testset$lat == uniqueLocations$lat[mindist.ind] & testset$lon == uniqueLocations$lon[mindist.ind])
    modRFI <- rep(0,length(smcts.ind))
    sevRFI <- rep(0,length(smcts.ind))
    modRFI[which((as.vector(sapply(testset$aquSFLAG[smcts.ind],function(x){ as.integer(intToBits(x))})))[29] == 1 |
                 (as.vector(sapply(testset$aquSFLAG[smcts.ind],function(x){ as.integer(intToBits(x))})))[31]) == 1] <- 1
    sevRFI[which((as.vector(sapply(testset$aquSFLAG[smcts.ind],function(x){ as.integer(intToBits(x))})))[30] == 1 |
                 (as.vector(sapply(testset$aquSFLAG[smcts.ind],function(x){ as.integer(intToBits(x))})))[32]) == 1] <- 1
    smcts <- data.frame(time=testset$time[smcts.ind], 
                        lat=testset$lat[smcts.ind],
                        lon=testset$lon[smcts.ind],
                        smc=testset$smc[smcts.ind],
                        modRFI=modRFI,
                        sevRFI=sevRFI)
#     locations_frame$name[station_ind] <- as.character(val_points$station[station_ind])
#     locations_frame$st_lat[station_ind] <- val_points$lat[station_ind]
#     locations_frame$st_lon[station_ind] <- val_points$lon[station_ind]
#     locations_frame$aqu_lat[station_ind] <- uniqueLocations$lat[mindist.ind]
#     locations_frame$aqu_lon[station_ind] <- uniqueLocations$lon[mindist.ind]
    #print(as.character(val_points$station[station_ind]))
  } else {
    smcts <- data.frame(time=-1,lat=-1,lon=-1,smc=-1)
  }
  
  filename = paste("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/timeseries/noRFImask2_final/era/",
                  val_points$station[station_ind], ".txt", sep="")
  write.table(smcts, 
             file=filename,
             quote=F,
             sep=",",
             row.names=F)
  
}

#write.table(locations_frame, file=locations_filename, quote=F, sep=",", row.names=F)
