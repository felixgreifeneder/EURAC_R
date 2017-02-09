# this scrip reads data for argentina and performs SVR training

library(raster)
library(gbm)
library(caret)
library(e1071)
library(ggplot2)

#read table and average measurements from same day at same location
features <- read.csv("X:/Workspaces/GrF/Processing/RS2_paper_argentina/RS2_argentina/buffered_features.csv", stringsAsFactors = F)
#trunc time from acquisition date
features$FECHA_HORA <- as.character(as.Date(features$FECHA_HORA, format = "%d/%m/%Y"))
#convert to matrix
features_mat <- as.matrix(features, dimnames=names(features))

uniqs <- unique(features_mat[,c(3,4)])
features_avg <- matrix(data = NA, nrow=nrow(uniqs), ncol=18)

#average redundant measurements
for (irow in 1:nrow(uniqs)){
  
  tmp_subset <- features_mat[which(features_mat[,3] == uniqs[irow,1] & features_mat[,4] == uniqs[irow,2]),]
  if (is.matrix(tmp_subset) == T){
    mean_smc <- as.character(mean(as.numeric(tmp_subset[,5])))
    features_avg[irow,] <- tmp_subset[1,]
    features_avg[irow,5] <- mean_smc
  } else
    features_avg[irow,] <- tmp_subset
  
}

#convert to dataframe and convert relevant rows to numeric
features_avg <- as.data.frame(features_avg, stringsAsFactors=F)
names(features_avg) <- names(features)
features_avg$RSOILMOIST <- as.numeric(features_avg$RSOILMOIST)
features_avg$hv <- as.numeric(features_avg$hv)
features_avg$stats_hv_m <- as.numeric(features_avg$stats_hv_m)
features_avg$stats_hv_s <- as.numeric(features_avg$stats_hv_s)
features_avg$hh <- as.numeric(features_avg$hh)
features_avg$stats_hh_m <- as.numeric(features_avg$stats_hh_m)
features_avg$stats_hh_s <- as.numeric(features_avg$stats_hh_s)
features_avg$slp <- as.numeric(features_avg$slp)
features_avg$asp <- as.numeric(features_avg$asp)
features_avg$hgt <- as.numeric(features_avg$hgt)
features_avg$ndvi <- as.numeric(features_avg$ndvi)
features_avg$incmean <- as.numeric(features_avg$incmean)
features_avg <- features_avg[,c(5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17)]

#convert to dB
features_avg$hv <- 10*log10(features_avg$hv)
features_avg$hh <- 10*log10(features_avg$hh)
#3x3 mean
features_avg$stats_hv_m <- 10*log10(features_avg$stats_hv_m)
features_avg$stats_hh_m <- 10*log10(features_avg$stats_hh_m)
#standard deviation
features_avg$stats_hv_s <- 10*log10(features_avg$stats_hv_s)
features_avg$stats_hh_s <- 10*log10(features_avg$stats_hh_s)

#filter feature set
val <- which(is.finite(features_avg$RSOILMOIST) & is.finite(features_avg$hv) &
               is.finite(features_avg$slp) & is.finite(features_avg$asp) & 
               is.finite(features_avg$hgt) & is.finite(features_avg$ndvi) & 
               is.finite(features_avg$hh) & is.finite(features_avg$incmean))
features_avg <- features_avg[val,]


#--------------------------------------------------------------------------------------------------
#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

set.seed(Sys.time())
trainIndices <- createDataPartition(features_avg$RSOILMOIST, p=0.8, list=FALSE, times=1)
training <- features_avg[trainIndices,]
testing <- features_avg[-trainIndices,]
#training <- CombDs
#testing <- CombDs

#--------------------------------------------------------------------------------------------------
#Training
#--------------------------------------------------------------------------------------------------
#SVR: package e1071

#setting the tuning options
#SVRtuningsettings <- tune.control(sampling = "cross", 
#                                  sampling.aggregate = mean,
#                                  sampling.dispersion = sd, 
#                                  cross = 5,
#                                  best.model = TRUE,
#                                  performances = TRUE)

#tuning the model (finding the best set  hyper parameters)
#SVRtuning <- tune(svm, RSOILMOIST ~ stats_hh_m + stats_hv_m + hgt + slp + asp + ndvi,
#                  data=training, 
#                  kernel="radial",
#                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
#                  tunecontrol=SVRtuningsettings)
#
#tunedModel <- SVRtuning$best.model

# save(SVRtuning, file="./Radarsat2Mazia/svrModel.dat")


#CARET
ctrl <- trainControl(method="cv", 
                     number = 10, 
                     search = "grid",
                     preProcOptions = c('center','scale'),
                     allowParallel = T)
mod <- train(RSOILMOIST ~ stats_hh_m + stats_hv_m + incmean + ndvi + hgt + slp + asp,
             data=training,
             method='svmRadial',
             tuneLength=5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)


#------------------------------------------------------------------------------------------------------------
#TESTING
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(mod, testing)
#SMCpredicted <- parallel_predictions(tunedModel, testing)
error <- sqrt(mean((testing$RSOILMOIST - SMCpredicted)^2))
r2 <- cor(testing$RSOILMOIST, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

tmp <- data.frame(x=testing$RSOILMOIST, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=4, height=3, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,0.5) + ylim(0,0.5) + coord_fixed(ratio = 1) +
      xlab('True SMC [m3m-3] \n') + ylab('\nEstimated SMC [m3m-3]') + 
      annotate("text", x=0.4, y=0.05, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(sqrt(r2), digits=3)), size=4)
p1
rm(tmp)