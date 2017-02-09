# this scrip reads data for argentina and performs SVR training

library(raster)
library(gbm)
library(caret)
library(e1071)
library(ggplot2)

# ARGENTINA

#read table and average measurements from same day at same location
features <- read.csv("C:/Users/FGreifeneder/Documents/tmp_proc/RS2_argentina/buffered_features.csv", stringsAsFactors = F)
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

features_argentina <- features_avg

remove(features_avg, features, features_mat)

#TERENO

#read table and average measurements from same day at same location
features_tereno <- read.csv("C:/Users/FGreifeneder/Documents/tmp_proc/RS2_Tereno/merged_dataset.csv", stringsAsFactors = F)


#filter feature set
val <- which(features_tereno$WC != 0)
features_tereno$WC <- features_tereno$WC/100


#Resize datasets so they have the same number of samples
Indices <- createDataPartition(features_argentina$RSOILMOIST, p=((nrow(features_argentina)-nrow(features_tereno))/nrow(features_argentina)), list=FALSE, times=1)
features_argentina <- features_argentina[-Indices,]


#merge datasets

features <- data.frame(smc=c(features_argentina$RSOILMOIST, features_tereno$WC),
                       hh=c(features_argentina$hh, features_tereno$hh),
                       hv=c(features_argentina$hv, features_tereno$hv),
                       lia=c(features_argentina$incmean, features_tereno$inc),
                       hgt=c(features_argentina$hgt, features_tereno$hgt),
                       asp=c(features_argentina$asp, features_tereno$asp),
                       slp=c(features_argentina$slp, features_tereno$slp),
                       ndvi=c(features_argentina$ndvi, features_tereno$ndvi))


#--------------------------------------------------------------------------------------------------
#partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
#be selected as training data. See - [http://topepo.github.io/caret/splitting.html]
print("Partitioning data in training and independent test set ...")

set.seed(Sys.time())
trainIndices <- createDataPartition(features$smc, p=0.8, list=FALSE, times=1)
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
mod <- train(smc ~ hh + hv + lia + ndvi + hgt + slp + asp,
             data=training,
             method='svmRadial',
             tuneLength = 5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=10), sigma= 10^seq(-2,2,len=10)),
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