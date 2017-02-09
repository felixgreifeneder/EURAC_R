#extract values for validation
library("e1071")
library("ggplot2")
library("R.matlab")
library("zoo")
source("./Aquarius/parallel_predictions.R")
library("raster")

#load required data
#--------------------------------------------------------------------------------------------------
val_points <- read.table("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_coordlist.txt", 
                         header=TRUE, sep=",")

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/train_test_mean_slope.dat")
load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/ascat_aqu/SVRmodel_performance.dat")
globcover <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/GLOBCOVER_L4_200901_200912_V2.3.tif")

#find insitu files
insitu.filelist = list.files(path = "C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_mariette",
                             pattern="*.csv",
                             full.names=T)
outdir = "C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_aqu_fltd/"


#predict
#--------------------------------------------------------------------------------------------------
tunedModel <- SVRtuning$best.model
 
SMCpredicted <- parallel_predictions(tunedModel, testset)


#find the matching time series
#--------------------------------------------------------------------------------------------------
# locations_filename <- "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/timeseries/locations.txt"
# locations_frame <- data.frame(name=character(nrow(val_points)), st_lat=numeric(nrow(val_points)), st_lon=numeric(nrow(val_points)), aqu_lat=numeric(nrow(val_points)), aqu_lon=numeric(nrow(val_points)))

merged_tss <- list()
allData <- data.frame()

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
    
    df.era <- data.frame(time=testset$time[smcts.ind], 
                          lat=testset$lat[smcts.ind],
                          lon=testset$lon[smcts.ind],
                          erasmc=testset$smc[smcts.ind])
    
    
    #initiate time series
    smcts.insitu <- read.zoo(insitu.filelist[station_ind], header=T, sep=",", 
                             format="%Y-%m-%d", aggregate=T, tz="GMT", index.column = 1)
    index(smcts.insitu) <- as.Date(index(smcts.insitu))
    smcts.pred <- zoo(SMCpredicted[smcts.ind], order.by=as.POSIXct(df.era$time*24*60*60, 
                                                                   origin="1900-01-01 00:00:00", 
                                                                   tz="GMT"))
    index(smcts.pred) <- as.Date(index(smcts.pred))
    smcts.era <- zoo(df.era$erasmc, order.by=as.POSIXct(df.era$time*24*60*60, 
                                                        origin="1900-01-01 00:00:00", 
                                                        tz="GMT"))
    index(smcts.era) <- as.Date(index(smcts.era))
    if ((length(unique(index(smcts.era))) != length(index(smcts.era))) |
      (length(unique(index(smcts.pred))) != length(index(smcts.pred))) |
      (length(unique(index(smcts.insitu$mean))) != length(index(smcts.insitu$mean)))){
      smcts.merged <- -1
      break
    }
  
    #merge time-series
    smcts.merged <- merge(smcts.pred, smcts.era, smcts.insitu$mean, all=c(T, T, F))
    
    #add to all-data array
    tmp <- fortify(smcts.merged)
    tmp <- data.frame(tmp, 
                      lc=rep(as.numeric(extract(globcover, 
                                                matrix(c(uniqueLocations$lon[mindist.ind], uniqueLocations$lat[mindist.ind]), 1, 2),
                                                buffer=1000,
                                                fun=modal,
                                                na.rm=T)),
                             nrow(tmp)))
    tmp <- data.frame(tmp, bias = tmp$smcts.pred-tmp$smcts.insitu.mean)
    allData <- rbind(allData, tmp)
    
    
  } else {
    smcts.merged <- -1
  }
  
  merged_tss[[station_ind]] <- smcts.merged
  
}

save(merged_tss, allData, file=paste(outdir,"tss.dat", sep=""))


#plot correlation
diff = allData$smcts.era - allData$smcts.insitu.mean
v_era <- which(abs(diff) < 0.05 & allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210)
allData <- allData[v_era,]

allData$lc <- as.factor(allData$lc)
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=smcts.pred, y=smcts.insitu.mean, colour=lc), data=allData) + 
  geom_point() + 
  #facet_wrap(~lc, ncol=3) +
  xlab("\nEstimated SMC [m3m-3]") + ylab("In-Situ SMC [m3m-3] \n") + 
  scale_colour_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Closed to open mixed forest",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Bare areas")) +
  theme(aspect.ratio=1) 
ggsave(filename=paste(outdir,"truevsest.png", sep=""))

#plot era vs insitu
allData$lc <- as.factor(allData$lc)
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=smcts.era, y=smcts.insitu.mean, colour=lc), data=allData[which(allData$lc != 220 & allData$lc != 210),]) + 
#ggplot(aes(x=smcts.era, y=smcts.insitu.mean, colour=lc), data=allData[v_era,]) + 
  geom_point() + 
  xlab("\nERA-Land SMC [m3m-3]") + ylab("In-Situ SMC [m3m-3] \n") + 
  scale_colour_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Closed brodleaved deciduous forest",
                                 "Closed needleleaved forest",
                                 "Open needleleaved forest",
                                 "Closed to open mixed forest",
                                 "Grassland (50-70%)/ forest or shrubland (20-50%)",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Sparse vegetation",
                                 "Bare areas")) +
 theme(aspect.ratio=1) 
ggsave(filename=paste(outdir,"truevsest.png", sep=""))

v_era <- which(abs(diff) < 0.05 & allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210)
allData <- allData[v_era,]
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=smcts.era, y=smcts.insitu.mean, colour=lc), data=allData[v_era,]) + 
  geom_point() + 
  xlab("\nERA-Land SMC [m3m-3]") + ylab("In-Situ SMC [m3m-3] \n") + 
  scale_colour_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Open needleleaved forest",
                                 "Closed to open mixed forest",
                                 "Grassland (50-70%)/ forest or shrubland (20-50%)",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Sparse vegetation",
                                 "Bare areas")) +
  theme(aspect.ratio=1) 
ggsave(filename=paste(outdir,"truevsest_filter.png", sep=""))

#plot bias
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=bias, fill=lc), data=allData) + 
  geom_boxplot() + xlab("\nLand-Cover") + ylab("Bias [m3m-3]\n") + ylim(-0.2,0.4) +
  scale_fill_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Closed to open mixed forest",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Bare areas"))
ggsave(filename=paste(outdir,"boxplot_bias.png", sep=""))

# plot era vs predicted

dev.new(width=7, height=5, noRStudioGD=T)
# ggplot(aes(x=smcts.era, y=smcts.insitu.mean, colour=lc), data=allData[which(allData$lc != 220 & allData$lc != 210),]) + 
ggplot(aes(x=smcts.era, y=smcts.insitu.mean, colour=lc), data=allData[v_era,]) + 
  geom_point() + 
  xlab("\nERA-Land SMC [m3m-3]") + ylab("In-Situ SMC [m3m-3] \n") + 
  scale_colour_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Closed brodleaved deciduous forest",
                                 "Closed needleleaved forest",
                                 "Open needleleaved forest",
                                 "Closed to open mixed forest",
                                 "Grassland (50-70%)/ forest or shrubland (20-50%)",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Sparse vegetation",
                                 "Bare areas")) +
  theme(aspect.ratio=1) 
ggsave(filename=paste(outdir,"truevsest.png", sep=""))

v_era <- which(allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210 & allData$lc != 200)
#allData <- allData[v_era,]
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=smcts.era, y=smcts.pred, colour=lc), data=allData[v_era,]) + 
  geom_point() + 
  xlab("\nERA-Land SMC [m3m-3]") + ylab("In-Situ SMC [m3m-3] \n") + 
  scale_colour_discrete(name="Land-Cover\n", 
                        labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                 "Vegetation (50-70%)/cropland (50-70%)",
                                 "Open needleleaved forest",
                                 "Closed to open mixed forest",
                                 "Grassland (50-70%)/ forest or shrubland (20-50%)",
                                 "Closed to open shrubland",
                                 "Closed to open grassland",
                                 "Sparse vegetation")) +
  theme(aspect.ratio=1) 
ggsave(filename=paste(outdir,"ERAvsest_filter.png", sep=""))


#metrics insitu vs pred
sink(paste(outdir, "metrics.txt", sep=""))

cat("Overall correlation")
cat("\n")
cat(cor(allData$smcts.insitu.mean, allData$smcts.pred, use="pairwise.complete.obs"))
cat("\n")
cat("Correlations by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(cor(allData$smcts.insitu.mean[userows], allData$smcts.pred[userows], use="pairwise.complete.obs"))
  cat("\n")
}
cat("\n")
cat("Overall mean bias")
cat("\n")
cat(mean(allData$smcts.pred-allData$smcts.insitu.mean, na.rm=T))
cat("\n")
cat("Mean Bias by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(mean(allData$smcts.pred[userows]-allData$smcts.insitu.mean[userows], na.rm=T))
  cat("\n")
}
cat("\n")
cat("Overall RMSE")
cat("\n")
cat(sqrt(mean((allData$smcts.pred-allData$smcts.insitu.mean)^2, na.rm=T)))
cat("\n")
cat("RMSE by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(sqrt(mean((allData$smcts.pred[userows]-allData$smcts.insitu.mean[userows])^2, na.rm=T)))
  cat("\n")
}
cat("\n")



sink()

#metrics era vs pred
sink(paste(outdir, "metrics_era.txt", sep=""))

cat("Overall correlation")
cat("\n")
cat(cor(allData$smcts.era, allData$smcts.pred, use="pairwise.complete.obs"))
cat("\n")
cat("Correlations by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(cor(allData$smcts.era[userows], allData$smcts.pred[userows], use="pairwise.complete.obs"))
  cat("\n")
}
cat("\n")
cat("Overall mean bias")
cat("\n")
cat(mean(allData$smcts.pred-allData$smcts.era, na.rm=T))
cat("\n")
cat("Mean Bias by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(mean(allData$smcts.pred[userows]-allData$smcts.era[userows], na.rm=T))
  cat("\n")
}
cat("\n")
cat("Overall RMSE")
cat("\n")
cat(sqrt(mean((allData$smcts.pred-allData$smcts.era)^2, na.rm=T)))
cat("\n")
cat("RMSE by LC")
cat("\n")
for (i in levels(allData$lc)){
  userows <- which(allData$lc == i)
  cat(i)
  cat(": ")
  cat(sqrt(mean((allData$smcts.pred[userows]-allData$smcts.era[userows])^2, na.rm=T)))
  cat("\n")
}
cat("\n")



sink()

# #plot time series
# for (i in 1:length(merged_tss)){
#   names(merged_tss[[i]]) <- c('Estimated', 'ERA-Land', 'In-situ')
#   ts_plot <- ggplot(aes(x= Index, y = Value, group = Series, colour = Series), data=fortify(merged_tss[[i]], melt= T)) +
#              geom_line() + xlab("Date") + ylab("") + facet_grid(Series ~ .) + theme(legend.position = "none") + geom_point()
#   
#   ggsave(filename=paste("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/merged_ts/Plot",i, ".png", sep=""), plot=ts_plot)
# }
