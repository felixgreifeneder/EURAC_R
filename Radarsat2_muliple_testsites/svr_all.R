# this scrip reads data for argentina and performs SVR training

library(raster)
library(gbm)
library(caret)
library(e1071)
library(ggplot2)

#read table and average measurements from same day at same location
features <- read.csv("C:/Users/FGreifeneder/Documents/tmp_proc/RS2_argentina_tereno/merged_dataset.csv", stringsAsFactors = F)
load("./Radarsat2_muliple_testsites/ComDs.dat")
features_Mazia <- CombDs


#filter feature set (Tereno and Argentina)
val <- which(features$WC != 0)
features <- features[val,]

#filter feature set (Mazia)
val <- which(features_Mazia$smc > 0 & features_Mazia$smc < 50 & features_Mazia$lia > 10 & features_Mazia$lia < 50)
features_Mazia <-  features_Mazia[val,]
features_Mazia$smc <- features_Mazia$smc / 100

#combine features sets
features_all <- data.frame(smc = c(features_Mazia$smc, features$WC),
                           hh = c(features_Mazia$hh, features$hh),
                           hv = c(features_Mazia$hv, features$hv),
                           lia = c(features_Mazia$lia, features$inc),
                           hgt = c(features_Mazia$height, features$hgt),
                           asp = c(features_Mazia$aspect, features$asp),
                           slp = c(features_Mazia$slope, features$slp),
                           ndvi = c(features_Mazia$ndvi, features$ndvi))


#--------------------------------------------------------------------------------------------------
#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

set.seed(Sys.time())
trainIndices <- createDataPartition(features_all$smc, p=0.8, list=FALSE, times=1)
training <- features_all#[trainIndices,]
testing <- features_all#[-trainIndices,]
#training <- CombDs
#testing <- CombDs

#--------------------------------------------------------------------------------------------------
#Training
#--------------------------------------------------------------------------------------------------
#CARET
ctrl <- trainControl(method="cv", 
                     number = 10, 
                     search = "grid",
                     preProcOptions = c('center','scale'),
                     allowParallel = T)
mod <- train(smc ~ hh + hv + lia + ndvi + hgt + slp + asp,
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
error <- sqrt(mean((testing$smc - SMCpredicted)^2))
r2 <- cor(testing$smc, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=testing$smc, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,0.5) + ylim(0,0.5) + coord_fixed(ratio = 1) + 
      xlab('True SMC [m3m-3] \n') + ylab('\nEstimated SMC [m3m-3]') + 
      annotate("text", x=0.4, y=0.05, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(sqrt(r2), digits=3)), size=4)
p1
rm(tmp)