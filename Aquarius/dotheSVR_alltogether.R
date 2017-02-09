#--------------------------------------------------------------------------------------------------
#loading required packages
library("e1071")
library("caret")
library("ggplot2")
library("reshape2")
library("hexbin")
library("RColorBrewer")
library("gridExtra")
source("./Aquarius/parallel_predictions.R")
library("raster")
library("Metrics")


#------------------------------------------------------------------------------------------------------------
#TUNING
#------------------------------------------------------------------------------------------------------------


#loading merged Aquarius/ASCAT/SMC dataset
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius/parameters_global_R")
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/train_test_mean_slope.dat")



#--------------------------------------------------------------------------------------------------
#visualizing the training set
#density plot training vs full data set
dev.new(width=8, height=6, noRStudioGD=T)
forPloting <- list(full=dataCollection$smc, training=trainset$smc)
long = melt(forPloting)
names(long) <- c("data", "Set")
ggplot(long, aes(x=data, colour=Set)) + geom_density(size=1) + ggtitle("SMC: full dataset vs. trainset") + 
  theme(plot.title=element_text(size=14, vjust=2, face="bold"))

#correlation plot training vs full data set
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=10, height=6, noRStudioGD=T)
p1 <- hexbinplot(smc~ascatM, data=dataCollection, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="Full dataset: SMC vs ascat mid-beam")
p2 <- hexbinplot(smc~ascatM, data=trainset, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="Training set: SMC vs ascat mid-beam")

grid.arrange(p1,p2,ncol=2)

#--------------------------------------------------------------------------------------------------
#Building the SVR model
#package e1071

#setting the tuning options
SVRtuningsettings <- tune.control(sampling = "cross", 
                                  sampling.aggregate = mean,
                                  sampling.dispersion = sd, 
                                  cross = 5,
                                  best.model = TRUE,
                                  performances = TRUE)

#tuning the model (finding the best set  hyper parameters)
SVRtuning <- tune(svm, smc ~ ascatF + ascatFmean + ascatFsd + 
                             ascatM + ascatMmean + ascatMsd +
                             ascatA + ascatAmean + ascatAsd +
                             aquVV + aquVVmean + aquVVsd +
                             aquVH + aquVHmean + aquVHsd,
                  data=trainset, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel <- SVRtuning$best.model


#caret

ctrl <- trainControl(method="cv", 
                     number = 5, 
                     search = "grid",
                     preProcOptions = c('center','scale'),
                     allowParallel = F)
mod <- train(smc ~ ascatF + ascatFmean + ascatFsd + 
               ascatM + ascatMmean + ascatMsd +
               ascatA + ascatAmean + ascatAsd,
             data=trainset,
             method='svmRadial',
             tuneLength=5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)

tunedModel <- mod

#------------------------------------------------------------------------------------------------------------
#TESTING
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
#SMCpredicted <- predict(tunedModel, testset)
SMCpredicted <- parallel_predictions(tunedModel, testset)
error <- rmse(testset$smc, SMCpredicted)
#error <- sqrt(mean((testset$smc - SMCpredicted)^2))
r <- cor(testset$smc, y=SMCpredicted, method="pearson")
print(paste("Error:",error,"R:",r))

tmp <- data.frame(x=testset$smc, y=SMCpredicted)
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- hexbinplot(y~x, data=tmp, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="True vs. Estimated SMC", xlab="true", ylab="estimated")
p1
rm(tmp)

#--------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------
#smc prediction for each point on the grid
#---------------------------------------------------------------------------------------------------------------

uniqueLocations <- unique(testset[,1:2])
nUL <- nrow(uniqueLocations)

prediction_performance <- data.frame(lat=rep(0,nUL), lon=rep(0,nUL), 
                                     rmse=rep(0,nUL),
                                     mse=rep(0,nUL),
                                     vari=rep(0,nUL),
                                     r2=rep(0,nUL), 
                                     bias=rep(0,nUL),
                                     n=rep(0,nUL), lc=rep(0,nUL))

globcover <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/GLOBCOVER_L4_200901_200912_V2.3.tif")

for (LocInd in c(1:nUL)){
  
  subset_indices <- which(testset$lat==uniqueLocations$lat[LocInd] & testset$lon==uniqueLocations$lon[LocInd])
  test_subset <- data.frame(testset[subset_indices,])
  
  SMCpredicted_subset <- SMCpredicted[subset_indices]
  #SMCpredicted <- predict(tunedModel, test_subset[,c("ascatM", "ascatincM", "aquVV", "aquVH", "lc")])
  #SMCpredicted <- SMC_predicted_full[subset_indices]
  rmse <- sqrt(mean((test_subset$smc - SMCpredicted_subset)^2))
  mse <- mean((test_subset$smc - SMCpredicted_subset)^2)
  vari <- var((test_subset$smc - SMCpredicted_subset)^2)
  r2 <- (cor(test_subset$smc, y=SMCpredicted_subset, method="pearson"))^2
  bias <- mean(test_subset$smc - SMCpredicted_subset)
  n <- length(subset_indices)
  
  prediction_performance$lat[LocInd] <- uniqueLocations$lat[LocInd]
  prediction_performance$lon[LocInd] <- uniqueLocations$lon[LocInd]
  prediction_performance$rmse[LocInd] <- rmse
  prediction_performance$mse[LocInd] <- mse
  prediction_performance$vari[LocInd] <- vari
  prediction_performance$r2[LocInd] <- r2
  prediction_performance$bias[LocInd] <- bias
  prediction_performance$n[LocInd] <- n
  prediction_performance$lc[LocInd] <- as.numeric(extract(globcover, matrix(c(uniqueLocations$lon[LocInd], uniqueLocations$lat[LocInd]), 1, 2)))
  
  #prediction_performance$lc[LocInd] <- test_subset$lc[1]
#   prediction_performance$mean_slope[LocInd] <- mean(test_subset$slope)
#   prediction_performance$mean_VHVV[LocInd] <- mean((10^(test_subset$aquVH/10)/10^(test_subset$aquVV/10)))
  
  
}

#save results
save(SVRtuning, prediction_performance, file="C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/ascat_aqu/SVRmodel_performance.dat")



#plot global grid with accuracies
map <- borders("world", colour="gray50", fill="gray50")
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
           geom_point(aes(x=lon, y=lat, colour=rmse), data=prediction_performance, size=1.2) + 
           scale_colour_gradientn(colours=r, name="RMSE", limit=c(0,0.3), na.value="red", breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("SVR prediction - RMSE\n")

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=r2), data=prediction_performance, size=1.2) + scale_colour_gradientn(colours=rev(r), name="R2")
gridPlt + ggtitle("SVR prediction - R2\n")

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=bias), data=prediction_performance, size=1.2) + 
          scale_colour_gradientn(colours=rev(r), name="Bias", limit=c(-0.4,0.4), breaks=c(-0.4,-0.2,0,0.2,0.4), labels=c("-0.4","-0.2","0","0.2","0.4"))
gridPlt + ggtitle("SVR prediction - Bias\n")


#boxplots-----------------------------------------------------

prediction_performance$lc <- as.factor(prediction_performance$lc)
new_lc <- prediction_performance$lc
new_lc <- as.character(new_lc)
new_lc[which(prediction_performance$lc == 14 | prediction_performance$lc == 20 | prediction_performance$lc == 30)] <- "croplands"
new_lc[which(prediction_performance$lc == 40 | prediction_performance$lc == 60 | prediction_performance$lc == 90 | prediction_performance$lc == 100)] <- "forest, open"
new_lc[which(prediction_performance$lc == 50 | prediction_performance$lc == 70)] <- "forest, closed"
new_lc[which(prediction_performance$lc == 110 | prediction_performance$lc == 120 | prediction_performance$lc == 130 | prediction_performance$lc == 140)] <- "shrublands/grasslands"
new_lc[which(prediction_performance$lc == 150 | prediction_performance$lc == 200)] <- "sparsely vegetated/bare"
new_lc[which(prediction_performance$lc == 160 | prediction_performance$lc == 170 | prediction_performance$lc == 180)] <- "regularly flooded"
new_lc[which(prediction_performance$lc == 190)] <- "urban"

prediction_performance <- data.frame(prediction_performance, lc_labels=new_lc)

use_lc <- which((prediction_performance$lc == 14 | prediction_performance$lc == 20 | prediction_performance$lc == 30 | prediction_performance$lc == 40 | 
                   prediction_performance$lc == 60 | prediction_performance$lc == 90 | prediction_performance$lc == 100 | prediction_performance$lc == 50 | 
                   prediction_performance$lc == 70 | prediction_performance$lc == 110 | prediction_performance$lc == 120 | prediction_performance$lc == 130 | 
                   prediction_performance$lc == 140 | prediction_performance$lc == 150 | prediction_performance$lc == 160 | prediction_performance$lc == 170 | 
                   prediction_performance$lc == 200 | prediction_performance$lc == 190 | prediction_performance$lc == 180) & 
                  prediction_performance$lat != -9999)
val <- use_lc

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(prediction_performance[val,], aes(x=lc_labels, y=r2, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(0,1) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("R2") + 
  ggtitle("R2\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(prediction_performance[val,], aes(x=lc_labels, y=rmse, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(0,0.4) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("RMSE") + 
  ggtitle("RMSE\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))


#package caret
# library(doSNOW)
# cl<-makeCluster(4, type="SOCK")
# registerDoSNOW(cl)
# 
# source("Aquarius/svr_define.R")
# 
# fitControl <- trainControl(method="cv",
#                            number=5,
#                            p=0.75,
#                            classProbs=FALSE,
#                            allowParallel=TRUE,
#                            verboseIter=TRUE,
#                            savePredictions=TRUE)
# 
# grdDef <- expand.grid(sigma=10^seq(-2,1,len=3), C=10^seq(-2,2,len=3), eps=10^seq(-2,-0.5,len=3), degree=1:2)
# 
# 
# SVRfitascat <- train(smc ~ ascatM + ascatincM,
#                      data=trainset,
#                      method=SVManova,
#                      trControl = fitControl,
#                      preProc=c("center", "scale"),
#                      tuneLength=20,
#                      allowParallel=TRUE,
#                      verbose=TRUE,
#                      tuneGrid=grdDef)
# 
# SMCpredicted <- predict.train(SVRfit, trainset[,c("aquVV", "aquVH", "aquINC")])
# error <- sqrt(mean((trainset$smc - SMCpredicted)^2))
# r2 <- cor(trainset$smc, y=SMCpredicted, method="pearson")





