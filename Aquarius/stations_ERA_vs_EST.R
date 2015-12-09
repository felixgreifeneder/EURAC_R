#extract values for validation
library("e1071")
library("ggplot2")
library("R.matlab")

#load required data
#--------------------------------------------------------------------------------------------------
val_points <- read.table("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/datapoints_SCA_validation.csv", 
                         header=TRUE, sep=",")
load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/testset_noRFIfilterR")
load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/SVR_ascatM_aqu/SVR_R")
testset <- testset_noRFIfilter
rm(testset_noRFIfilter)

ascatM_quant <- quantile(testset$ascatM)
ascatM_IQR <- IQR(testset$ascatM)
aquVV_quant <- quantile(testset$aquVV)
aquVV_IQR <- IQR(testset$aquVV)
aquVH_quant <- quantile(testset$aquVH)
aquVH_IQR <- IQR(testset$aquVH)

valid <- which(testset$ascatM > ascatM_quant[2] - 1.5*ascatM_IQR &
                 testset$ascatM < ascatM_quant[4] + 1.5*ascatM_IQR &
                 testset$aquVV > aquVV_quant[2] - 1.5*aquVV_IQR &
                 testset$aquVV < aquVV_quant[4] + 1.5*aquVV_IQR &
                 testset$aquVH > aquVH_quant[2] - 1.5*aquVH_IQR &
                 testset$aquVH < aquVH_quant[4] + 1.5*aquVH_IQR)

testset <- testset[valid,]

# Cdata <- readMat("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/SCA_Eumetsat_Claudia/new test data noRFI/estim_smc_vv_f_vhvv_ts_n2.mat")
# testset <- as.data.frame(Cdata$datats.n[,2:17])
# colnames(testset) <- c("lat",
#                        "lon",
#                        "time",
#                        "ascatF",
#                        "ascatincF",
#                        "ascatM",
#                        "ascatincM",
#                        "ascatA",
#                        "ascatincA",
#                        "slope",
#                        "aquVV",
#                        "aquVH",
#                        "aquINC",
#                        "aquSFLAG",
#                        "smc",
#                        "TEMP")
# 
# SMCpredicted <- t(Cdata$estim.smc.vv.f.vhvv.ts.n2)

#plot
# map <- borders("world", colour="gray50", fill="gray50")
# dev.new(width=10, height=6, noRStudioGD=T)
# gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=n), data=prediction_performance, size=1.2)
# gridPlt <- gridPlt + geom_point(aes(x=lon, y=lat), data=val_points, colour="red", size=2)
# gridPlt + ggtitle("Global grid\n")


#predict
#--------------------------------------------------------------------------------------------------
tunedModel <- SVRtuning$best.model

SMCpredicted <- predict(tunedModel, testset)


#find the matching time series
#--------------------------------------------------------------------------------------------------
smcts <- data.frame(time=numeric(0), lat=numeric(0), lon=numeric(0), smc=numeric(0), era=numeric(0), modRFI=numeric(0), sevRFI=numeric(0), st=numeric(0))

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
    smcts <- rbind(smcts, data.frame(time=testset$time[smcts.ind], 
                                     lat=testset$lat[smcts.ind],
                                     lon=testset$lon[smcts.ind],
                                     smc=SMCpredicted[smcts.ind],
                                     era=testset$smc[smcts.ind],
                                     modRFI=modRFI,
                                     sevRFI=sevRFI,
                                     st=station_ind))
  } 
}

smcts$st <- as.factor(smcts$st)
valid <- which(smcts$sevRFI != 1 & smcts$modRFI != 1)
ggplot(smcts[valid,], aes(x=smc, y=era, color=RFI)) + geom_point(shape=1, size=3) + scale_x_continuous(limits=c(0,0.5)) + scale_y_continuous(limits=c(0,0.5))
map <- borders("world", colour="gray50", fill="gray50")
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=st), data=smcts[valid,], size=3)
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt
