#this routine applies an iterative training approach selecting a good subset of each training and re-training on the rest

source("Aquarius/readAquInfo.R")
library("scales")
library("ggplot2")
library("e1071")
library("caret")
library("hexbin")
library("RColorBrewer")
library("gridExtra")
library("matlab")

load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/train_testset_pseudRand_R")


# trainset_use <- which((10^(trainset$aquVH/10)/10^(trainset$aquVV/10)) > 0.1 & (10^(trainset$aquVH/10)/10^(trainset$aquVV/10)) < 0.4)
# trainset <- trainset[trainset_use, ]
# testset_use <- which((10^(testset$aquVH/10)/10^(testset$aquVV/10)) > 0.1 & (10^(testset$aquVH/10)/10^(testset$aquVV/10)) < 0.4)
# testset <- testset[testset_use, ]

newTrainset <- data.frame(lat=numeric(), 
                          lon=numeric(), 
                          smc=numeric(), 
                          ascatM=numeric(), 
                          ascatincM=numeric(), 
                          aquVV=numeric(), 
                          aquVH=numeric(),
                          trainIT = numeric())

trainingIter <- 0

repeat {

  #-----------------------------------------------------SVR-------------------------------------------------------
  trainingIter <- trainingIter+1
  print(trainingIter)
  
  #package e1071
  SVRtuningsettings <- tune.control(sampling = "cross", 
                                    sampling.aggregate = mean,
                                    sampling.dispersion = sd, 
                                    cross = 5,
                                    best.model = TRUE,
                                    performances = TRUE)
  
  #scale the training set
  # trainset_scaled <- data.frame(smc=(trainset$smc-min(trainset$smc))/(max(trainset$smc)-min(trainset$smc)),
  #                               ascatM=(trainset$ascatM-min(trainset$ascatM))/(max(trainset$ascatM)-min(trainset$ascatM)))
  
  SVRtuning <- tune(svm, smc ~ ascatM + ascatincM + aquVV + aquVH, 
                    data=trainset, 
                    kernel="radial",
                    ranges = list(epsilon = logspace(-2,-0.5,n=3), gamma=logspace(-2,1,n=3), cost=logspace(-2,2,n=3)),
                    #ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-4,4,len=10)),
                    tunecontrol=SVRtuningsettings, verbose=TRUE)
  

  tunedModel <- SVRtuning$best.model
  
  print("Overall performance based on training set:")
  SMCpredicted <- predict(tunedModel, trainset[,c("ascatM", "ascatincM", "aquVV", "aquVH")])
  error <- sqrt(mean((trainset$smc - SMCpredicted)^2))
  r2 <- cor(trainset$smc, y=SMCpredicted, method="pearson")^2
  print(paste("Error:",error,"R2:",r2))
  
  #claculate performance for each grid point
  uniqueLocations <- unique(trainset[,1:2])
  nUL <- nrow(uniqueLocations)
  
#   prediction_performance <- data.frame(lat=rep(0,nUL), lon=rep(0,nUL), rmse=rep(0,nUL), r2=rep(0,nUL))
#   
#   for (LocInd in c(1:nUL)){
#     
#     subset_indices <- which(trainset$lat==uniqueLocations$lat[LocInd] & trainset$lon==uniqueLocations$lon[LocInd])
#     train_subset <- data.frame(trainset[subset_indices,])
#     
#     #discard potins with too little data
#     if (length(train_subset$smc) <= 50) {
#       r2 <- -1
#       error <- -1
#     } else {
#       SMCpredicted <- predict(tunedModel, train_subset)
#       #SMCpredicted <- SMC_predicted_full[subset_indices]
#       error <- sqrt(mean((train_subset$smc - SMCpredicted)^2))
#       r2 <- (cor(train_subset$smc, y=SMCpredicted, method="pearson"))^2
#     }
#     
#     prediction_performance$lat[LocInd] <- uniqueLocations$lat[LocInd]
#     prediction_performance$lon[LocInd] <- uniqueLocations$lon[LocInd]
#     prediction_performance$rmse[LocInd] <- error
#     prediction_performance$r2[LocInd] <- r2
#     #   prediction_performance$mean_slope[LocInd] <- mean(test_subset$slope)
#     #   prediction_performance$mean_VHVV[LocInd] <- mean((10^(test_subset$aquVH/10)/10^(test_subset$aquVV/10)))
#     
#     
#   }
  
  #filter locations with rmse < 0.05, the rest is used for re-training
#   good_locations <- which(prediction_performance$rmse < 0.05 & prediction_performance$r2 > 0.4)
#   good_samples <- numeric()
#   if (length(good_locations) > 0){
#     for (LocInd in good_locations){
#       good_samples <- c(good_samples, which(trainset$lat==prediction_performance$lat[LocInd] & trainset$lon==prediction_performance$lon[LocInd]))
#     }
#     
#     if (length(good_samples >= 50)){
#       newTrainset <- rbind(newTrainset, data.frame(trainset[good_samples,c("lat","lon","smc","ascatM","ascatincM","aquVV","aquVH")], trainIT=rep(trainingIter, length(good_samples))))
#       trainset <- trainset[-good_samples,]  
#     }
#   }

  tmp <- abs(trainset$smc - SMCpredicted)
  good_samples <- which(tmp < 0.05)
  if (length(good_samples) >= 50){
    newTrainset <- rbind(newTrainset, data.frame(trainset[good_samples,c("lat","lon","smc","ascatM","ascatincM","aquVV","aquVH")], trainIT=rep(trainingIter, length(good_samples))))
    trainset <- trainset[-good_samples,]   
  } else {
    break
  }
  
  #if (nrow(trainset) < 200 | length(good_samples) < 50) break

}

#--------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------

#define the classification training set. 
classTrainset <- rbind(newTrainset, data.frame(trainset[,c("lat","lon","smc","ascatM","ascatincM","aquVV","aquVH")], trainIT=rep(trainingIter, nrow(trainset))))
classTrainset$trainIT <- factor(classTrainset$trainIT)


#re-train SVR, using iteration labels from first training phase as an additional input


SVRtuning <- tune(svm, smc ~ ascatM + ascatincM + aquVV + aquVH + trainIT, 
                  data=newTrainset, 
                  kernel="radial",
                  ranges = list(epsilon = logspace(-2,-0.5,n=3), gamma=logspace(-2,1,n=3), cost=logspace(-2,2,n=3)),
                  #ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-4,4,len=10)),
                  tunecontrol=SVRtuningsettings)

tunedModel <- SVRtuning$best.model


print("Overall performance based on full testset:")
SMCpredicted <- predict(tunedModel, newTrainset[,c("ascatM", "ascatincM", "aquVV", "aquVH","trainIT")])
error <- sqrt(mean((newTrainset$smc - SMCpredicted)^2))
r2 <- cor(newTrainset$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=newTrainset$smc, y=SMCpredicted)
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- hexbinplot(y~x, data=tmp, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="True vs. Estimated SMC", xlab="true", ylab="estimated")
p1
rm(tmp)


#--------------------------------------------------------------------------------------------------------------
# Classification of test-set
#--------------------------------------------------------------------------------------------------------------

#apply iteration labels from training to the test set, based on nearest neighbour
#claculate performance for each grid point
uniqueTestLocations <- unique(testset[,1:2])
uniqueTrainingLocations <- unique(classTrainset[,1:2])
nUL <- nrow(uniqueTestLocations)
nUL_train <- nrow(uniqueTrainingLocations)

testset_trainIT <- rep(5, times=nrow(testset))

for (LocInd in c(1:nUL)){
  
  subset_indices <- which(testset$lat==uniqueTestLocations$lat[LocInd] & testset$lon==uniqueTestLocations$lon[LocInd])
  #train_subset <- data.frame(trainset[subset_indices,])
  
  distances <- numeric()
  for (LocInd_train in c(1:nUL_train)){
    distances <- c(distances, sqrt((uniqueTestLocations$lat[LocInd]-uniqueTrainingLocations$lat[LocInd_train])^2 + (uniqueTestLocations$lon[LocInd]-uniqueTrainingLocations$lon[LocInd_train])^2))
  }
  
  nn <- uniqueTrainingLocations[which.min(distances),]
  subset_training_indices <- which(classTrainset$lat==nn$lat & classTrainset$lon==nn$lon)
  testset_trainIT[subset_indices] <- classTrainset$trainIT[subset_training_indices[1]]
  
}

testset <- data.frame(testset, trainIT=factor(testset_trainIT))

#save results
save(newTrainset, classTrainset, testset, file="C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/SVR_iterative/testing_training_new_R")

#---------------------------------------------------------------------------------------------------------------
#SMC predictions
#---------------------------------------------------------------------------------------------------------------
#prediction on full testset

test_subset <- testset[which(testset$trainIT != 2),c("lat","lon","smc","ascatM","ascatincM","aquVV","aquVH","trainIT")]
test_subset$trainIT <- as.numeric(test_subset$trainIT)

print("Overall performance based on full testset:")
SMCpredicted <- predict(tunedModel, test_subset[,c("ascatM", "ascatincM", "aquVV", "aquVH","trainIT")])
error <- sqrt(mean((test_subset$smc - SMCpredicted)^2))
r2 <- cor(test_subset$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=test_subset$smc, y=SMCpredicted)
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- hexbinplot(y~x, data=tmp, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="True vs. Estimated SMC", xlab="true", ylab="estimated")
p1
rm(tmp)

#performance for each grid point
uniqueLocations <- unique(test_subset[,1:2])
nUL <- nrow(uniqueLocations)

prediction_performance <- data.frame(lat=rep(0,nUL), lon=rep(0,nUL), rmse=rep(0,nUL), r2=rep(0,nUL))

for (LocInd in c(1:nUL)){
  
  subset_indices <- which(test_subset$lat==uniqueLocations$lat[LocInd] & test_subset$lon==uniqueLocations$lon[LocInd])
  test_subset2 <- data.frame(test_subset[subset_indices,])
  
  SMCpredicted <- predict(tunedModel, test_subset2)
  #SMCpredicted <- SMC_predicted_full[subset_indices]
  error <- sqrt(mean((test_subset2$smc - SMCpredicted)^2))
  r2 <- (cor(test_subset2$smc, y=SMCpredicted, method="pearson"))^2
  
  prediction_performance$lat[LocInd] <- uniqueLocations$lat[LocInd]
  prediction_performance$lon[LocInd] <- uniqueLocations$lon[LocInd]
  prediction_performance$rmse[LocInd] <- error
  prediction_performance$r2[LocInd] <- r2
  #   prediction_performance$mean_slope[LocInd] <- mean(test_subset$slope)
  #   prediction_performance$mean_VHVV[LocInd] <- mean((10^(test_subset$aquVH/10)/10^(test_subset$aquVV/10)))
  
  
}

#save results
save(SVRtuning, prediction_performance, file="C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessingCRMASK/SVR_ascatM_inc_Aqu_OPTI_RMSE_R2/SVR_R")

#plot global grid with accuracies
map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
           geom_point(aes(x=lon, y=lat, colour=rmse), data=prediction_performance, size=1.2) + 
           scale_colour_gradientn(colours=r, name="RMSE", limit=c(0,0.3), na.value="red", breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("SVR prediction - RMSE\n")

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
           geom_point(aes(x=lon, y=lat, colour=r2), data=prediction_performance, size=1.2) + 
           scale_colour_gradientn(colours=rev(r), name="R2", limit=c(0,1))
gridPlt + ggtitle("SVR prediction - R2\n")
      





