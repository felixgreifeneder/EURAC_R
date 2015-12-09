#this routine computes correlations between smc and vv or vh, respectively

source("Aquarius/readAquInfo.R")
library("scales")
library("ggplot2")
library("e1071")
library("caret")
library("hexbin")
library("RColorBrewer")
library("gridExtra")
library("matlab")

load("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessing20150409/train_testset_pseudRand_R")


trainset_use <- which((10^(trainset$aquVH/10)/10^(trainset$aquVV/10)) > 0.1 & (10^(trainset$aquVH/10)/10^(trainset$aquVV/10)) < 0.4)
trainset <- trainset[trainset_use, ]
testset_use <- which((10^(testset$aquVH/10)/10^(testset$aquVV/10)) > 0.1 & (10^(testset$aquVH/10)/10^(testset$aquVV/10)) < 0.4)
testset <- testset[testset_use, ]


#-----------------------------------------------------SVR-------------------------------------------------------
#package e1071
SVRtuningsettings <- tune.control(sampling = "cross", 
                                  sampling.aggregate = mean,
                                  sampling.dispersion = sd, 
                                  cross = 5,
                                  best.model = TRUE,
                                  performances = TRUE)

#scale the training set
# trainset_scaled <- data.frame(smc=(trainset$smc-min(trainset$smc))/(max(trainset$smc)-min(trainset$smc)),
#                               ascatM=(trainset$ascatM-min(trainset$ascatM))/(max(trainset$ascatM)-min(trainset$ascatM)))

SVRtuning <- tune(svm, ascatM ~ ascatincM + smc, 
                  data=trainset, 
                  kernel="radial",
                  ranges = list(epsilon = logspace(-2,-0.5,n=3), gamma=logspace(-2,1,n=3), cost=logspace(-2,2,n=3)),
                  #ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-4,4,len=10)),
                  tunecontrol=SVRtuningsettings, verbose=TRUE)

#scale the testset
# testset_scaled <- data.frame(smc=(testset$smc-min(trainset$smc))/(max(trainset$smc)-min(trainset$smc)),
#                              ascatM=(testset$ascatM-min(trainset$ascatM))/(max(trainset$ascatM)-min(trainset$ascatM)))

tunedModel <- SVRtuning$best.model

write.svm(tunedModel, svm.file="C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessingCRMASK/SVR_ascatMinc_SMC/fittedModel.svm",
          scale.file="C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessingCRMASK/SVR_ascatMinc_SMC/scalingX.scale",
          yscale.file="C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessingCRMASK/SVR_ascatMinc_SMC/scalingY.yscale")

print("Overall performance based on trainingset:")
SIG0predicted <- predict(tunedModel, trainset[,c("ascatincM", "smc")])
error <- sqrt(mean((trainset$ascatM - SIG0predicted)^2))
r2 <- cor(trainset$ascatM, y=SIG0predicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

writetxtData(SIG0Predicted_out, "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/Aquarius/global_v2/beam2_reprocessingCRMASK/SVR_ascatMinc_SMC/predicted_sig0_trainingset.txt")

tmp <- data.frame(x=trainset$ascatM, y=SIG0predicted)
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
dev.new(width=8, height=6, noRStudioGD=T)
p1 <- hexbinplot(y~x, data=tmp, colramp=rf, 
                 aspect=1, type="r", trans=log, 
                 inv=exp, main="True vs. Estimated SIG0", xlab="true", ylab="estimated")
p1
rm(tmp)





