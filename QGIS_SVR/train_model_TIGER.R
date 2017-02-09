library(raster)
library(gbm)
library(caret)
library(ggplot2)

samples_file <- "X:\\Workspaces\\GrF\\Processing\\ESA_TIGER\\all_field_data.txt"
SIG0_VV <- T
PLIA <- F
S2_b2 <- F
S2_b3 <- F
S2_b4 <- F
S2_b5 <- F
S2_b6 <- F
S2_b7 <- F
S2_b8 <- F
Land_Cover <- T
Training <- 0.8
Training_Testing_Split <- "X:\\Workspaces\\GrF\\Processing\\ESA_TIGER\\QGIS_SVR_traintest.dat"
Output_Directory <- "X:\\Workspaces\\GrF\\Processing\\ESA_TIGER"


# Load samples from csv-file
samples <- read.csv(samples_file)

# construct samples data frame
samples.df <- data.frame(id=c(1:length(samples$Ave_Moist)), stringsAsFactors = F)
samples.df$Ave_Moist <- as.numeric(samples$Ave_Moist)
samples.df$Sig0_VV <- as.numeric(10*log10(samples$Sig0_VV_lin))
samples.df$PLIA <- as.numeric(samples$PLIA)
samples.df$S_Band_2 <- as.numeric(samples$S_Band_2)
samples.df$S_Band_3 <- as.numeric(samples$S_Band_3)
samples.df$S_Band_4 <- as.numeric(samples$S_Band_4)
samples.df$S_Band_5 <- as.numeric(samples$S_Band_5)
samples.df$S_Band_6 <- as.numeric(samples$S_Band_6)
samples.df$S_Band_7 <- as.numeric(samples$S_Band_7)
samples.df$S_Band_8 <- as.numeric(samples$S_Band_8)
samples.df$LC <- as.numeric(samples$LC)

# filter available samples
# check for NA values
validNumb <- rep(TRUE, nrow(samples.df))
for (i in c(1:ncol(samples.df))){
  
  validNumb <- validNumb & is.finite(samples.df[,i])
  
}
# filter LIA values < 10° and > 50°
validLia <- samples.df$PLIA >= 10 & samples.df$PLIA <= 50

# combined valid samples
valid <- validNumb & validLia

# check for invalid SM values
valid.sm <- samples.df$Ave_Moist > 0 & samples.df$Ave_Moist < 70
# combine
valid <- valid & valid.sm

samples.df <- samples.df[valid,]

# generate formula
SMC.forumla <- "Ave_Moist ~ "

if (SIG0_VV == T){
  SMC.forumla <- paste(SMC.forumla, 'Sig0_VV', sep='')
}
if (PLIA == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'PLIA', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + PLIA', sep='')
  }
}
if (S2_b2 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_2', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_2', sep='')
  }
}
if (S2_b3 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_3', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_3', sep='')
  }
}
if (S2_b4 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_4', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_4', sep='')
  }
}
if (S2_b5 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_5', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_5', sep='')
  }
}
if (S2_b6 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_6', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_6', sep='')
  }
}
if (S2_b7 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_7', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_7', sep='')
  }
}
if (S2_b8 == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'S_Band_8', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + S_Band_8', sep='')
  }
}
if (Land_Cover == T){
  if (SMC.forumla =="Ave_Moist ~ "){
    SMC.forumla <- paste(SMC.forumla, 'LC', sep='')
  } else {
    SMC.forumla <- paste(SMC.forumla, ' + LC', sep='')
  }
}


if (Training == 1){
  training <- samples.df
  testing <- samples.df
} else {
  if (samples_file != "") {
    #partitioning the data in training and testing. Due to the amoung of data available only a small subset shall
    print("Partitioning data in training and independent test set ...")
    
    startSet <- sample(1:nrow(samples.df), 5)
    samplePool <- as.matrix(samples.df[-startSet, 3:ncol(samples.df)])
    start <- as.matrix(samples.df[startSet, 3:ncol(samples.df)])
    newSamp <- maxDissim(start, samplePool, n=floor(nrow(samplePool)*Training))
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
mod <- train(as.formula(SMC.forumla), data=samples.df,
             method='svmRadial',
             tuneLength=5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)

save(mod, file = paste(Output_Directory,'/QGIS_SVR_model.dat', sep=''))


#Testing of trained model on the independent test-set
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(mod, testing)
#SMCpredicted <- parallel_predictions(tunedModel, testing)
error <- sqrt(mean((testing$Ave_Moist - SMCpredicted)^2))
r <- cor(testing$Ave_Moist, y=SMCpredicted, method="pearson")
print(paste("Error:",error,"R:",r))

tmp <- data.frame(x=testing$Ave_Moist, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=9, height=7, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,50) + ylim(0,50) + coord_fixed(ratio = 1) +
  xlab('True SMC [m3m-3] \n') + ylab('\nEstimated SMC [m3m-3]') + 
  annotate("text", x=40, y=10, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(r, digits=3)), size=4)
p1
rm(tmp)
ggsave(filename = paste(Output_Directory,'/True_vs_Estimated.png', sep=""), p1)


