# this scrip reads data for argentina and performs SVR training

library(raster)
library(gbm)
library(caret)
library(e1071)
library(ggplot2)

#read table and average measurements from same day at same location
features <- read.csv("X:/Workspaces/GrF/Processing/RS2_paper_argentina/RS2_Tereno/data20150716.csv", stringsAsFactors = F)


#filter feature set
val <- which(features$WC != 0)
features$WC <- features$WC/100


#--------------------------------------------------------------------------------------------------
#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

set.seed(Sys.time())
trainIndices <- createDataPartition(features$WC, p=0.8, list=FALSE, times=1)
training <- features#[trainIndices,]
testing <- features#[-trainIndices,]
#training <- CombDs
#testing <- CombDs

#--------------------------------------------------------------------------------------------------
#Training
#--------------------------------------------------------------------------------------------------
#CARET
ctrl <- trainControl(method="LOOCV", 
                             number = 10, 
                             search = "grid",
                             preProcOptions = c('center','scale'),
                             allowParallel = T)
mod <- train(WC ~ hh + hv + inc + ndvi + hgt + slp + asp,
             data=training,
             method='svmRadial',
             tuneLength = 5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)


#------------------------------------------------------------------------------------------------------------
#TESTING
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(mod, testing)
#SMCpredicted <- parallel_predictions(tunedModel, testing)
error <- sqrt(mean((testing$WC - SMCpredicted)^2))
r2 <- cor(testing$WC, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=testing$WC, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=4, height=3, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,0.50) + ylim(0,0.50) + coord_fixed(ratio = 1) +
      xlab('True SMC [m3m-3] \n') + ylab('\nEstimated SMC [m3m-3]') + 
      annotate("text", x=0.4, y=0.05, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(sqrt(r2), digits=3)), size=4)
p1
rm(tmp)