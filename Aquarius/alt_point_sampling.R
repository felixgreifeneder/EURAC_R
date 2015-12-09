source("Aquarius/readAquInfo.R")
source("Aquarius/write_data_to_textfile.R")
library("scales")
library("ggplot2")
library("e1071")
library("caret")
library("raster")

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius/parameters_global_R")

testset <- data.frame(time = numeric(),
                      lat = numeric(),
                      lon = numeric(),
                      ascatF = numeric(),
                      ascatincF = numeric(),
                      ascatFmean = numeric(),
                      ascatFsd = numeric(),
                      ascatM = numeric(),
                      ascatincM = numeric(),
                      ascatMmean = numeric(),
                      ascatMsd = numeric(),
                      ascatA = numeric(),
                      ascatincA = numeric(),
                      ascatAmean = numeric(),
                      ascatAsd = numeric(),
                      slope = numeric(),
                      slopemean = numeric(),
                      aquVV = numeric(),
                      aquVVmean = numeric(),
                      aquVVsd = numeric(),
                      aquVH = numeric(),
                      aquVHmean = numeric(),
                      aquVHsd = numeric(),
                      aquHH = numeric(),                      
                      aquHHmean = numeric(),
                      aquHHsd = numeric(),
                      aquINC = numeric(),
                      aquRFIV = numeric(),
                      aquRFIH = numeric(),
                      eraSMC = numeric(),
                      TEMP = numeric())

trainset <- testset


#filtering
print("Partitioning data in training and independent test set ...")

ascatM_quant <- quantile(dataCollection$ascatM)
ascatM_IQR <- IQR(dataCollection$ascatM)
aquVV_quant <- quantile(dataCollection$aquVV)
aquVV_IQR <- IQR(dataCollection$aquVV)
aquVH_quant <- quantile(dataCollection$aquVH)
aquVH_IQR <- IQR(dataCollection$aquVH)
aquHH_quant <- quantile(dataCollection$aquHH)
aquHH_IQR <- IQR(dataCollection$aquHH)

valid <- which(dataCollection$ascatM > ascatM_quant[2] - 1.5*ascatM_IQR &
               dataCollection$ascatM < ascatM_quant[4] + 1.5*ascatM_IQR &
               dataCollection$aquVV > aquVV_quant[2] - 1.5*aquVV_IQR &
               dataCollection$aquVV < aquVV_quant[4] + 1.5*aquVV_IQR &
               dataCollection$aquVH > aquVH_quant[2] - 1.5*aquVH_IQR &
               dataCollection$aquVH < aquVH_quant[4] + 1.5*aquVH_IQR &
               dataCollection$aquHH > aquHH_quant[2] - 1.5*aquHH_IQR &
               dataCollection$aquHH < aquHH_quant[4] + 1.5*aquVH_IQR)

dataCollection <- dataCollection[valid,]

#add rows for mean data
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))
dataCollection <- cbind(dataCollection, rep.int(0, nrow(dataCollection)))

names(dataCollection) <- c("lat","lon","time","ascatF","ascatincF","ascatM","ascatincM","ascatA","ascatincA",
                           "slope","aquVV","aquVH","aquHH","aquINC","aquSFLAG","smc","TEMP","aquVVmean","aquVVsd",
                           "aquVHmean","aquVHsd","aquHHmean","aquHHsd",
                           "ascatFmean", "ascatFsd",
                           "ascatMmean", "ascatMsd",
                           "ascatAmean", "ascatAsd", 
                           "slopemean", "slopesd")

#calculate the mean for sig0 for each grid point
uniqueLoc = unique(dataCollection[,c("lat", "lon")])
for (i in 1:nrow(uniqueLoc)){
  
  locIds <- which(dataCollection$lat == uniqueLoc[i,1] & dataCollection$lon == uniqueLoc[i,2])
  
  dataCollection$aquVVmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$aquVV[locIds]/10))), length(locIds))
  dataCollection$aquVHmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$aquVH[locIds]/10))), length(locIds)) 
  dataCollection$aquHHmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$aquHH[locIds]/10))), length(locIds))
  dataCollection$ascatFmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$ascatF[locIds]/10))), length(locIds))
  dataCollection$ascatMmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$ascatM[locIds]/10))), length(locIds))
  dataCollection$ascatAmean[locIds] <- rep.int(10*log10(mean(10^(dataCollection$ascatA[locIds]/10))), length(locIds))
  dataCollection$slopemean[locIds] <- rep.int(mean(dataCollection$slope[locIds]), length(locIds))
  if (length(locIds)>1){
    dataCollection$aquVVsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$aquVV[locIds]/10))), length(locIds))
    dataCollection$aquVHsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$aquVH[locIds]/10))), length(locIds))
    dataCollection$aquHHsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$aquHH[locIds]/10))), length(locIds))
    dataCollection$ascatFsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$ascatF[locIds]/10))), length(locIds))
    dataCollection$ascatMsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$ascatM[locIds]/10))), length(locIds))
    dataCollection$ascatAsd[locIds] <- rep.int(10*log10(sd(10^(dataCollection$ascatA[locIds]/10))), length(locIds))
    dataCollection$slopesd[locIds] <- rep.int(sd(dataCollection$slope[locIds]), length(locIds))
  } else {
    dataCollection$aquVVsd[locIds] <- 0
    dataCollection$aquVHsd[locIds] <- 0
    dataCollection$aquHHsd[locIds] <- 0
    dataCollection$ascatFsd[locIds] <- 0
    dataCollection$ascatMsd[locIds] <- 0
    dataCollection$ascatAsd[locIds] <- 0
    dataCollection$slopesd[locIds] <- 0
  }
  
  
  
}

dataCollection <- dataCollection[which(dataCollection$aquVVsd != 0),]


# 
# #plot global filtered grid with valid points
# uniqueLocations <- unique(filteredDataset[,1:2])
# uniqueIndices <- as.numeric(rownames(uniqueLocations))
map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat), data=dataCollection, colour="red", size=1.2)
gridPlt + ggtitle("Global grid")


#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

sizeofDataset <- dim(dataCollection)
sizeofDataset <- sizeofDataset[1]
#nTrain <- round(sizeofDataset*0.01)
nTrain <- 5000
pTrain <- 5000*(1/sizeofDataset)
#nTrain <- nTrain - 5
alltrainIndices <- NULL

#--------------------------------------------SAMPLING----------------------------------------------------

#create a new data frame only holding information for algorithm training (e.g.: no latlon) to select the subset
nfD <- data.frame(dataCollection$smc, 
                  dataCollection$ascatF,
                  dataCollection$ascatFmean,
                  dataCollection$ascatFsd,
                  dataCollection$ascatM, 
                  dataCollection$ascatMmean,
                  dataCollection$ascatMsd,
                  dataCollection$ascatA,
                  dataCollection$ascatAmean,
                  dataCollection$ascatAsd,
                  dataCollection$slope, 
                  dataCollection$slopemean,
                  dataCollection$slopesd,
                  dataCollection$aquVV, 
                  dataCollection$aquVVmean,
                  dataCollection$aquVVsd,
                  dataCollection$aquVH,
                  dataCollection$aquVHmean,
                  dataCollection$aquVHsd,
                  dataCollection$aquHH,
                  dataCollection$aquHHmean,
                  dataCollection$aquHHsd)
names(nfD) <- c("smc", 
                "ascatF",
                "ascatFmean",
                "ascatFsd", 
                "ascatM", 
                "ascatMmean",
                "ascatMsd",
                "ascatA",
                "ascatAmean",
                "ascatAsd",
                "slope", 
                "slopemean",
                "slopesd",
                "aquVV", 
                "aquVVmean", 
                "aquVVsd",
                "aquVH",
                "aquVHmean",
                "aquVHsd",
                "aquHH",
                "aquHHmean",
                "aquHHsd")

#pseudo random sampling (simple)
set.seed(3456)
trainIndices <- createDataPartition(nfD$smc, p=pTrain, list=FALSE, times=1)
trainset <- dataCollection[trainIndices,]
testset <- dataCollection[-trainIndices,]

  
#sampling based on max difference
# dataScaled <- scale(nfD)
# set.seed(5)
# startSet <- sample(c(1:nrow(dataScaled)), 15)
# samplePool <- dataScaled[-startSet,]
# start <- dataScaled[startSet,]
# newSamp <- maxDissim(start, samplePool, n=nTrain, randomFrac=.1, useNames=T, verbose=T)
#     
# trainIndices <- c(startSet, newSamp)
# 
# #testset <- rbind(testset, selectedRows[-trainIndices,])
# trainset <- dataCollection[trainIndices,]
# testset <- dataCollection[-trainIndices,]

save(trainset, testset, file="C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/train_test_mean_slope.dat")
#writetxtData(trainset, "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/training_txt.txt")
#writetxtData(testset, "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150417/testing_txt.txt")


#--------------------------------------------------------------------------------------------------------
#---------------------------------------------PLOTING----------------------------------------------------

#plot global grid based on the trainingset
uniqueLocations <- unique(trainset[,1:2])
uniqueIndices <- as.numeric(rownames(uniqueLocations))
map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat), data=trainset, colour="red", size=1.2)
gridPlt + ggtitle("Training Set")

