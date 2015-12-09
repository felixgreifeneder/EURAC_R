library("R.matlab")
library("hexbin")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")

Cdata <- readMat("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/SCA_Eumetsat_Claudia/new results/estim_smc_vv_f_slope_ts_n.mat")
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

estimates <- t(Cdata$estim.smc.vv.f.vhvv.ts.n)

print("Overall performance based on full testset:")
error <- sqrt(mean((testset$smc - estimates)^2))
r2 <- cor(testset$smc, y=estimates, method="pearson")^2
print(paste("Error:",error,"R2:",r2))


#get performance measures
#--------------------------------------------------------------------------------------------------


uniqueLocations <- unique(testset[,1:2])
nUL <- nrow(uniqueLocations)

prediction_performance <- data.frame(lat=rep(0,nUL), lon=rep(0,nUL), 
                                     rmse=rep(0,nUL),
                                     mse=rep(0,nUL),
                                     vari=rep(0,nUL),
                                     r2=rep(0,nUL), 
                                     n=rep(0,nUL))

for (LocInd in c(1:nUL)){
  
  subset_indices <- which(testset$lat==uniqueLocations$lat[LocInd] & testset$lon==uniqueLocations$lon[LocInd])
  test_subset <- data.frame(testset[subset_indices,])
  
  SMCpredicted <- estimates[subset_indices,]
  #SMCpredicted <- SMC_predicted_full[subset_indices]
  rmse <- sqrt(mean((test_subset$smc - SMCpredicted)^2))
  mse <- mean((test_subset$smc - SMCpredicted)^2)
  vari <- var((test_subset$smc - SMCpredicted)^2)
  r2 <- (cor(test_subset$smc, y=SMCpredicted, method="pearson"))^2
  n <- length(subset_indices)
  
  prediction_performance$lat[LocInd] <- uniqueLocations$lat[LocInd]
  prediction_performance$lon[LocInd] <- uniqueLocations$lon[LocInd]
  prediction_performance$rmse[LocInd] <- rmse
  prediction_performance$mse[LocInd] <- mse
  prediction_performance$vari[LocInd] <- vari
  prediction_performance$r2[LocInd] <- r2
  prediction_performance$n[LocInd] <- n
  #   prediction_performance$mean_slope[LocInd] <- mean(test_subset$slope)
  #   prediction_performance$mean_VHVV[LocInd] <- mean((10^(test_subset$aquVH/10)/10^(test_subset$aquVV/10)))
  
  
}


#plot global grid with accuracies
#---------------------------------------------------------------------------------------------------
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=rmse), data=prediction_performance, size=1.2) + 
  scale_colour_gradientn(colours=r, name="RMSE", limit=c(0,0.3), na.value="red", breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("RMSE: ASCAT+Aquarius\n")

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=r2), data=prediction_performance, size=1.2) + scale_colour_gradientn(colours=rev(r), name="R2")
gridPlt + ggtitle("R2: ASCAT+Aquarius\n")

tmp <- data.frame(x=testset$smc, y=estimates)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- hexbinplot(y~x, data=tmp, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="True vs. Estimated SMC", xlab="true", ylab="estimated")
p1
rm(tmp)