#attempt estimation with SVR using in-situ data for training

library("rhdf5")
library("ggplot2")
library("gtable")
library("grid")
library("e1071")
library("caret")
source("./SMAP/dailySMCavg.R")

#create load SIG0 time-series

load("./SMAP/SIG0TSs.dat")

B1SIG0$time <- as.Date(B1SIG0$time)
B2SIG0$time <- as.Date(B2SIG0$time)
B3SIG0$time <- as.Date(B3SIG0$time)


b1 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB10.dat")
b2 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB15.dat")
b3 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB20.dat")
# b1 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP10_scaled.dat")
# b2 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP15_scaled.dat")
# b3 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP20_scaled.dat")

b1$date <- as.Date(b1$date)
b2$date <- as.Date(b2$date)
b3$date <- as.Date(b3$date)

#--------------------------------------------------------------------------------------------------
#B1

SScols1 <- data.frame(date=as.character(),
                     smc=as.numeric(),
                     sig0vv_fore=as.numeric(),
                     sig0vv_aft=as.numeric(),
                     sig0hh_fore=as.numeric(),
                     sig0hh_aft=as.numeric(),
                     sig0xpol_fore=as.numeric(),
                     sig0xpol_aft=as.numeric(), stringsAsFactors = F)
cntr <- 1
for (i in 1:length(B1SIG0$time)){
  if (B1SIG0$dist[i] >= 1 | B1SIG0$vv_fore[i]==-9999) {next}
  if (any(b1$date == B1SIG0$time[i])==F) {next}
  if (B1SIG0$rflag[i] >= 8 & B1SIG0$rflag[i] < 16) {next}
  SScols1[cntr,1] <- as.character.Date(B1SIG0$time[i])
  SScols1[cntr,2] <- b1$mean_SM05[b1$date == B1SIG0$time[i]]
  SScols1[cntr,3] <- B1SIG0$vv_fore[i]
  SScols1[cntr,4] <- B1SIG0$vv_aft[i]
  SScols1[cntr,5] <- B1SIG0$hh_fore[i]
  SScols1[cntr,6] <- B1SIG0$hh_aft[i]
  SScols1[cntr,7] <- B1SIG0$xpol_fore[i]
  SScols1[cntr,8] <- B1SIG0$xpol_aft[i]
  cntr <- cntr+1
  
}

#--------------------------------------------------------------------------------------------------
#B2

SScols2 <- data.frame(date=as.character(),
                     smc=as.numeric(),
                     sig0vv_fore=as.numeric(),
                     sig0vv_aft=as.numeric(),
                     sig0hh_fore=as.numeric(),
                     sig0hh_aft=as.numeric(),
                     sig0xpol_fore=as.numeric(),
                     sig0xpol_aft=as.numeric(), stringsAsFactors = F)
cntr <- 1
for (i in 1:length(B2SIG0$time)){
  if (B2SIG0$dist[i] >= 1 | B2SIG0$vv_fore[i]==-9999) {next}
  if (any(b2$date == B2SIG0$time[i])==F) {next}
  if (B2SIG0$rflag[i] >= 8 & B2SIG0$rflag[i] < 16) {next}
  SScols2[cntr,1] <- as.character.Date(B2SIG0$time[i])
  SScols2[cntr,2] <- b2$mean_SM05[b2$date == B2SIG0$time[i]]
  SScols2[cntr,3] <- B2SIG0$vv_fore[i]
  SScols2[cntr,4] <- B2SIG0$vv_aft[i]
  SScols2[cntr,5] <- B2SIG0$hh_fore[i]
  SScols2[cntr,6] <- B2SIG0$hh_aft[i]
  SScols2[cntr,7] <- B2SIG0$xpol_fore[i]
  SScols2[cntr,8] <- B2SIG0$xpol_aft[i]
  cntr <- cntr+1
  
}

#--------------------------------------------------------------------------------------------------
#B3

SScols3 <- data.frame(date=as.character(),
                     smc=as.numeric(),
                     sig0vv_fore=as.numeric(),
                     sig0vv_aft=as.numeric(),
                     sig0hh_fore=as.numeric(),
                     sig0hh_aft=as.numeric(),
                     sig0xpol_fore=as.numeric(),
                     sig0xpol_aft=as.numeric(), stringsAsFactors = F)
cntr <- 1
for (i in 1:length(B3SIG0$time)){
  if (B3SIG0$dist[i] >= 1 | B3SIG0$vv_fore[i]==-9999) {next}
  if (any(b1$date == B3SIG0$time[i])==F) {next}
  if (B2SIG0$rflag[i] >= 8 & B2SIG0$rflag[i] < 16) {next}
  SScols3[cntr,1] <- as.character.Date(B3SIG0$time[i])
  SScols3[cntr,2] <- b3$mean_SM05[b3$date == B3SIG0$time[i]]
  SScols3[cntr,3] <- B3SIG0$vv_fore[i]
  SScols3[cntr,4] <- B3SIG0$vv_aft[i]
  SScols3[cntr,5] <- B3SIG0$hh_fore[i]
  SScols3[cntr,6] <- B3SIG0$hh_aft[i]
  SScols3[cntr,7] <- B3SIG0$xpol_fore[i]
  SScols3[cntr,8] <- B3SIG0$xpol_aft[i]
  cntr <- cntr+1
  
}

#--------------------------------------------------------------------------------------------------
#SVR estimation test
#reshaping sig0vsstation data frame

SVRdt <- rbind(SScols1, SScols2, SScols3)


#Splitting the dataset into training and test-set

nrofrows <- nrow(SVRdt)

#setting the size of the training set
nTrain <- 23
pTrain <- 23*(1/nrofrows)

#using a pseudo random sampling strategy
set.seed(3456)
trainIndices <- createDataPartition(SVRdt$smc, p=pTrain, list=FALSE, times=1)
#trainIndices <- trainIndices[,1]
#load("./SMAP/trainingindices.dat")
save(trainIndices, file="./SMAP/trainingindices.dat")
training <- SVRdt[trainIndices,]
testing <- SVRdt[-trainIndices,]

#--------------------------------------------------------------------------------------------------
#Building the SVR model
#package e1071

#setting the tuning options
SVRtuningsettings <- tune.control(sampling = "fix", 
                                  fix = 2/3,
                                  best.model = TRUE,
                                  performances = TRUE)

#tuning the model (finding the best set  hyper parameters)
SVRtuning <- tune(svm, smc ~ sig0vv_fore + sig0vv_aft + sig0hh_fore + sig0hh_aft + sig0xpol_fore + sig0xpol_aft, 
                  data=training, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol = SVRtuningsettings)

tunedModel <- SVRtuning$best.model

#--------------------------------------------------------------------------------------------------
#Estimating SMC

print("Overall performance based on full testset:")
SMCpredicted <- predict(tunedModel, testing)
error <- sqrt(mean((testing$smc - SMCpredicted)^2))
r2 <- cor(testing$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(true=testing$smc, predicted=SMCpredicted)

p2 <- ggplot(data=tmp, aes(x=true,y=predicted)) + 
  geom_point()  +
  theme(aspect.ratio=1) +
  geom_abline(intercept = 0, slope = 1, linetype=2) + 
  xlim(0,50) + ylim(0,50)

dev.new()
p2