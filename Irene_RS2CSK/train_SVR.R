library(raster)
library(gbm)
library(caret)
library(ggplot2)
library(R.matlab)

samples_file <- readMat("X:\\Workspaces\\NiI\\05_Code\\MOMS_new\\input\\mazia_data_past.mat")
samples.df <- as.data.frame(samples_file$matVect)
samples.df <- data.frame(id=c(1:nrow(samples.df)), 
                         V1=samples.df$V1,
                         V2=samples.df$V2,
                         V3=samples.df$V3,
                         V4=samples.df$V4,
                         V5=samples.df$V5,
                         V6=samples.df$V6,
                         V7=samples.df$V7,
                         V8=samples.df$V8)
Training <- 0.8
Training_Testing_Split <- ""
Output_Directory <- "X:\\Workspaces\\NiI\\05_Code\\MOMS_new\\output_felix"

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
    save(training_ids, testing_ids, file = paste(Output_Directory, '/traintest_mazia_data_past.dat', sep=''))
    
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
mod <- train(training[,3:ncol(training)], training[,2],
             method='svmRadial',
             tuneLength=5,
             #tuneGrid = data.frame(C=10^seq(-2,2,len=20), sigma= 10^seq(-2,1,len=20)),
             trControl=ctrl)

save(mod, file = paste(Output_Directory,'/SVR_model_mazia_data_past.dat', sep=''))


#Testing of trained model on the independent test-set
#------------------------------------------------------------------------------------------------------------



print("Overall performance based on full testset:")
SMCpredicted <- predict(mod, testing[,3:ncol(testing)])
#SMCpredicted <- parallel_predictions(tunedModel, testing)
bias <- mean(testing$V1 - SMCpredicted)
mse <- mean((testing$V1 - SMCpredicted)^2)
error <- sqrt(mean((testing$V1 - SMCpredicted)^2))
r <- cor(testing$V1, y=SMCpredicted, method="pearson")
print(paste("Error:",error,"R:",r))

#output results to testfile
sink(paste(Output_Directory, '/report_mazia_data_past.txt', sep=""))
cat("===============================================\n")
cat("SVR Training Report\n")
cat("===============================================\n")
cat("\n")
cat("Input: mazia_data_past.mat\n")
cat("Size of training set:", nrow(training), '\n')
cat("Size of test set:", nrow(testing), '\n')
cat("\n")
cat("Accuracies based on test set:\n")
cat("MSE=", mse, ",", "RMSE=", error, ",", "Bias=", bias, ",", "R=", r, "\n")
cat("------------------------------------------------\n")
sink()



# create plot

tmp <- data.frame(x=testing$V1, y=SMCpredicted)
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=9, height=7, noRStudioGD=T)
p1 <- ggplot(data=tmp, aes(x=x, y=y)) + geom_point() +geom_smooth(method=lm) + xlim(0,70) + ylim(0,70) + coord_fixed(ratio = 1) +
  xlab('True SMC [Vol. %] \n') + ylab('\nEstimated SMC [Vol. %]') + 
  annotate("text", x=60, y=10, label=paste('RMSE: ', format(error, digits=3), '\n', 'R: ', format(r, digits=3)), size=4)
p1
rm(tmp)
ggsave(filename = paste(Output_Directory,'/tr_vs_est_mazia_data_past.png', sep=""), p1)



