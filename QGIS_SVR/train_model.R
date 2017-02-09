library(raster)
library(gbm)
library(caret)
library(ggplot2)

samples_file <- "X:/Workspaces/GrF/Processing/SVR_QGIS/testfile_for_qgis.csv"
SIG0_VV <- F
SIG0_VH <- F
SIG0_HH <- T
SIG0_HV <- T
LIA <- F
Height <- F
Slope <- F
Aspect <- F
NDVI <- F
Land_Cover <- F
Training <- 0.8
Training_Testing_Split <- ""
Output_Directory <- "X:/Workspaces/GrF/Processing/SVR_QGIS"


# Load samples from csv-file
samples <- read.csv(samples_file)

# construct samples data frame
samples.df <- data.frame(id=c(1:length(samples$SMC)), stringsAsFactors = F)
samples.df$smc <- as.numeric(samples$SMC)
if (SIG0_VV == T){
  samples.df$VV <- as.numeric(samples$vv)
}
if (SIG0_VH == T){
  samples.df$VH <- as.numeric(samples$vh)
}
if (SIG0_HH == T){
  samples.df$HH <- as.numeric(samples$hh)
}
if (SIG0_HV == T){
  samples.df$HV <- as.numeric(samples$hv)
}
if (LIA == T){
  samples.df$LIA <- as.numeric(samples$lia)
}
if (Height == T){
  samples.df$HGT <- as.numeric(samples$hgt)
}
if (Slope == T){
  samples.df$SLP <- as.numeric(samples$slp)
}
if (Aspect == T){
  samples.df$ASP <- as.numeric(samples$asp)
}
if (NDVI == T){
  samples.df$NDVI <- as.numeric(samples$ndvi)
}
if (Land_Cover == T){
  samples.df$LC <- as.numeric(samples$lc)
}

# filter available samples
# check for NA values
validNumb <- rep(TRUE, nrow(samples.df))
for (i in c(1:ncol(samples.df))){
  
  validNumb <- validNumb & is.finite(samples.df[,i])
  
}
# filter LIA values < 10° and > 50°
validLia <- samples.df$LIA >= 10 & samples.df$LIA <= 50
# valid topography
validHgt <- samples.df$HGT > 0
validAsp <- samples.df$ASP >= 0 & samples.df$ASP <= 360
validSlp <- samples.df$SLP > 0

# combined valid samples
valid <- validNumb & validLia & validHgt & validAsp & validSlp
samples.df <- samples.df[valid,]
  
if (samples_file != "") {
  #partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
  print("Partitioning data in training and independent test set ...")
  
  startSet <- sample(1:nrow(samples.df), 5)
  samplePool <- as.matrix(samples.df[-startSet, 3:ncol(samples.df)])
  start <- as.matrix(samples.df[startSet, 3:ncol(samples.df)])
  newSamp <- maxDissim(start, samplePool, n=floor(nrow(samplePool)*0.8))
  # extract training and testing set
  samplePool <- samples.df[-startSet, 1:ncol(samples.df)]
  training <- rbind(samples.df[startSet, 1:ncol(samples.df)], samplePool[newSamp,])
  testing <- samplePool[-newSamp,]
  training_ids <- training$id
  testing_ids <- testing$id
  save(training_ids, testing_ids, file = paste(Output_Directory, '/QGIS_SVR_traintest.dat', sep=''))
  
} else {
  
  load(Training_Testing_Split)
  
  tind <- c()
  for (i in training_ids){
    tind <- c(tind, which(samples.df$id == i))
  }
  training <- samples.df[tind,]
  tind <- c()
  for (i in testing_ids){
    tind <- c(tind, which(samples.df$id == i))
  }
  testing <- samples.df[tind,]
  
  
}

# Train SVR model
# -----------------------------------------------------------------------------
set.seed(Sys.time())

# Definition of training parameters
ctrl <- trainControl(method="cv", 
                     number = 10, 
                     search = "grid",
                     preProcOptions = c('center','scale'),
                     allowParallel = T)

# Model training based on "training" dataset
mod <- train(training[,3:ncol(training)], training[,2],
             method='svmRadial',
             tuneLength=5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)

save(mod, file = paste(Output_Directory,'/QGIS_SVR_model.dat', sep=''))


#Testing of trained model on the independent test-set
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(mod, testing[,3:ncol(testing)])
#SMCpredicted <- parallel_predictions(tunedModel, testing)
error <- sqrt(mean((testing$smc - SMCpredicted)^2))
r <- cor(testing$smc, y=SMCpredicted, method="pearson")
print(paste("Error:",error,"R:",r))

tmp <- data.frame(x=testing$smc, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=9, height=7, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,0.5) + ylim(0,0.5) + coord_fixed(ratio = 1) +
  xlab('True SMC [m3m-3] \n') + ylab('\nEstimated SMC [m3m-3]') + 
  annotate("text", x=0.4, y=0.05, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(r, digits=3)), size=4)
p1
rm(tmp)
ggsave(filename = paste(Output_Directory,'/True_vs_Estimated.png', sep=""), p1)


