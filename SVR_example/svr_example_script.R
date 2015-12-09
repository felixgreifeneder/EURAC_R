#Initialisation
#--------------------------------------------------------------------------------------------------
#loading required packages
library("e1071")
library("caret")
library("ggplot2")
library("reshape2")
library("hexbin")
library("RColorBrewer")
library("gridExtra")


#loading merged Aquarius/ASCAT/SMC dataset
load("./SVR_example/dataset_australiaR")


#plot the data locations
#map <- borders("world", region="australia", colour="gray50", fill="gray50")
map <- get_map(location=c(111,-47,155,-8.5), zoom=4, maptype="terrain", source="google")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggmap(map) + 
           geom_point(aes(x=lon, y=lat), data=data_australia, colour="red", size=2)
gridPlt + ggtitle("Data grid\n")

#--------------------------------------------------------------------------------------------------
#Splitting the dataset into training and test-set

nrofrows <- nrow(data_australia)

#setting the size of the training set
nTrain <- 1000
pTrain <- 1000*(1/nrofrows)

#using a pseudo random sampling strategy
set.seed(3456)
trainIndices <- createDataPartition(data_australia$smc, p=pTrain, list=FALSE, times=1)
training <- data_australia[trainIndices,]
testing <- data_australia[-trainIndices,]

#--------------------------------------------------------------------------------------------------
#visualizing the training set
#density plot training vs full data set
dev.new(width=8, height=6, noRStudioGD=T)
forPloting <- list(full=data_australia$smc, training=training$smc)
long = melt(forPloting)
names(long) <- c("data", "Set")
ggplot(long, aes(x=data, colour=Set)) + geom_density(size=1) + ggtitle("SMC: full dataset vs. trainset") + 
  theme(plot.title=element_text(size=14, vjust=2, face="bold"))

#correlation plot training vs full data set
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=10, height=6, noRStudioGD=T)
p1 <- hexbinplot(smc~ascatM, data=data_australia, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="Full dataset: SMC vs ascat mid-beam")
p2 <- hexbinplot(smc~ascatM, data=training, colramp=rf, 
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
SVRtuning <- tune(svm, smc ~ ascatM + ascatincM + aquVV + aquVH, 
                  data=training, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel <- SVRtuning$best.model

#--------------------------------------------------------------------------------------------------
#Estimating SMC

print("Overall performance based on full testset:")
SMCpredicted <- predict(tunedModel, testing)
error <- sqrt(mean((testing$smc - SMCpredicted)^2))
r2 <- cor(testing$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

#Plot true vs estimated SMC
tmp <- data.frame(x=testing$smc, y=SMCpredicted)
dev.new(width=8, height=6, noRStudioGD=T)
truevsest <- hexbinplot(y~x, data=tmp, colramp=rf, 
                        aspect=1, type="r", trans=log, 
                        inv=exp, main="True vs. Estimated SMC", xlab="true", ylab="estimated")
truevsest
rm(tmp)

#--------------------------------------------------------------------------------------------------
#SMC estimation for each grid point individually

uniqueLocations <- unique(testing[,1:2])
nL <- nrow(uniqueLocations)

prediction_performance <- data.frame(lat=rep(0,nL), lon=rep(0,nL), 
                                     rmse=rep(0,nL),
                                     mse=rep(0,nL),
                                     vari=rep(0,nL),
                                     r2=rep(0,nL), 
                                     n=rep(0,nL))

for (LocInd in c(1:nL)){
  
  subset_indices <- which(testing$lat==uniqueLocations$lat[LocInd] & testing$lon==uniqueLocations$lon[LocInd])
  test_subset <- data.frame(testing[subset_indices,])
  
  SMCpredicted_subset <- SMCpredicted[subset_indices]
  #SMCpredicted <- SMC_predicted_full[subset_indices]
  rmse <- sqrt(mean((test_subset$smc - SMCpredicted_subset)^2))
  mse <- mean((test_subset$smc - SMCpredicted_subset)^2)
  vari <- var((test_subset$smc - SMCpredicted_subset)^2)
  r2 <- (cor(test_subset$smc, y=SMCpredicted_subset, method="pearson"))^2
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

#plotting the results
#plot global grid with accuracies
map <- get_map(location=c(111,-47,155,-8.5), zoom=4, maptype="terrain", source="google")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggmap(map) + 
           geom_point(aes(x=lon, y=lat, colour=rmse), data=prediction_performance, size=2) + 
           scale_colour_gradientn(colours=r, name="RMSE", 
                                  limit=c(0,0.3), na.value="red", 
                                  breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("SVR prediction - RMSE\n")

map <- get_map(location=c(111,-47,155,-8.5), zoom=4, maptype="terrain", source="google")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggmap(map) +
           geom_point(aes(x=lon, y=lat, colour=r2), data=prediction_performance, size=2) + 
           scale_colour_gradientn(colours=rev(r), name="R2", limit=c(0,1))
gridPlt + ggtitle("SVR prediction - R2\n")

