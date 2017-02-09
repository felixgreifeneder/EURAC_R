#train SVR model

library(gbm)
library(caret)
library(e1071)
library(ggplot2)
source("./Aquarius/parallel_predictions.R")

load("./Radarsat2Mazia/ComDs.dat")

#CombDs <- CombDs[,-8]

#filter dataset

val <- which(CombDs$smc > 0 & CombDs$smc < 50 & CombDs$lia > 10 & CombDs$lia < 50)
CombDs <- CombDs[val,]


#--------------------------------------------------------------------------------------------------
#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

set.seed(Sys.time())
trainIndices <- createDataPartition(CombDs$smc, p=0.8, list=FALSE, times=1)
training <- CombDs[trainIndices,]
testing <- CombDs[-trainIndices,]
#training <- CombDs
#testing <- CombDs

#--------------------------------------------------------------------------------------------------
#Training
#--------------------------------------------------------------------------------------------------
#SVR: package e1071

#setting the tuning options
SVRtuningsettings <- tune.control(sampling = "cross", 
                                  sampling.aggregate = mean,
                                  sampling.dispersion = sd, 
                                  cross = 5,
                                  best.model = TRUE,
                                  performances = TRUE)

#tuning the model (finding the best set  hyper parameters)
SVRtuning <- tune(svm, smc ~ hh + hv + lia + height + slope + aspect,
                  data=training, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel <- SVRtuning$best.model

save(SVRtuning, file="./Radarsat2Mazia/svrModel.dat")


#------------------------------------------------------------------------------------------------------------
#TESTING
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(tunedModel, testing)
#SMCpredicted <- parallel_predictions(tunedModel, testing)
error <- sqrt(mean((testing$smc - SMCpredicted)^2))
r2 <- cor(testing$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=testing$smc, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point()
p1
rm(tmp)