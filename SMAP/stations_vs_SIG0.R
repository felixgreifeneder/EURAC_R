#compare stations and radiometer BT

library("rhdf5")
library("ggplot2")
library("gtable")
library("grid")
library("e1071")
library("caret")
library("RColorBrewer")
source("./SMAP/dailySMCavg.R")

#create load SIG0 time-series

load("./SMAP/SIG0TSsInterpol.dat")

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


SIG0vsSMC <- data.frame(smc=as.numeric(), 
                        sig0f=as.numeric(), 
                        sig0a=as.numeric(),
                        pol=as.character(), 
                        stat=as.character(), 
                        date=as.character(), 
                        NDIST=as.numeric(), 
                        LIAf=as.numeric(),
                        LIAa=as.numeric(), stringsAsFactors = F)



cntr <- 1

for (i in 1:length(B1SIG0$time)){
  
  #if (B1SIG0$dist[i] >= 0.006) {next}
  if (any(b1$date == B1SIG0$time[i]) == F) {next}
  #if (B1SIG0$rflag[i] >= 8 & B1SIG0$rflag[i] < 16) {next}
  #print(i)
  SIG0vsSMC[cntr,1] <- b1$mean_SM05[which(b1$date == B1SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B1SIG0$vv_fore[i] + B1SIG0$vv_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B1SIG0$vv_fore[i]
  SIG0vsSMC[cntr,3] <- B1SIG0$vv_aft[i]
  SIG0vsSMC[cntr,4] <- "VV"
  SIG0vsSMC[cntr,5] <- "B1000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b1$date[which(b1$date == B1SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B1SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B1SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B1SIG0$lia_aft[i]
  
  cntr <- cntr + 1
  SIG0vsSMC[cntr,1] <- b1$mean_SM05[which(b1$date == B1SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B1SIG0$hh_fore[i] + B1SIG0$hh_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B1SIG0$hh_fore[i]
  SIG0vsSMC[cntr,3] <- B1SIG0$hh_aft[i]
  SIG0vsSMC[cntr,4] <- "HH"
  SIG0vsSMC[cntr,5] <- "B1000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b1$date[which(b1$date == B1SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B1SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B1SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B1SIG0$lia_aft[i]
  cntr <- cntr + 1

  SIG0vsSMC[cntr,1] <- b1$mean_SM05[which(b1$date == B1SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B1SIG0$xpol_fore[i] + B1SIG0$xpol_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B1SIG0$xpol_fore[i]
  SIG0vsSMC[cntr,3] <- B1SIG0$xpol_aft[i]
  SIG0vsSMC[cntr,4] <- "XPOL"
  SIG0vsSMC[cntr,5] <- "B1000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b1$date[which(b1$date == B1SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B1SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B1SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B1SIG0$lia_aft[i]
  cntr <- cntr + 1

  
}

for (i in 1:length(B2SIG0$time)){
  
  #if (B2SIG0$dist[i] >= 0.006) {next}
  if (any(b2$date == B2SIG0$time[i]) == F) {next}
  #if (B2SIG0$rflag[i] >= 8 & B2SIG0$rflag[i] < 16) {next}
  #print(i)
  SIG0vsSMC[cntr,1] <- b2$mean_SM05[which(b2$date == B2SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B2SIG0$vv_fore[i] + B2SIG0$vv_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B2SIG0$vv_fore[i]
  SIG0vsSMC[cntr,3] <- B2SIG0$vv_aft[i]
  SIG0vsSMC[cntr,4] <- "VV"
  SIG0vsSMC[cntr,5] <- "B1500"
  SIG0vsSMC[cntr,6] <- as.character.Date(b2$date[which(b2$date == B2SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B2SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B2SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B2SIG0$lia_aft[i]
  cntr <- cntr + 1

  SIG0vsSMC[cntr,1] <- b2$mean_SM05[which(b2$date == B2SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B2SIG0$hh_fore[i] + B2SIG0$hh_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B2SIG0$hh_fore[i]
  SIG0vsSMC[cntr,3] <- B2SIG0$hh_aft[i]
  SIG0vsSMC[cntr,4] <- "HH"
  SIG0vsSMC[cntr,5] <- "B1500"
  SIG0vsSMC[cntr,6] <- as.character.Date(b2$date[which(b2$date == B2SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B2SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B2SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B2SIG0$lia_aft[i]
  cntr <- cntr + 1

  SIG0vsSMC[cntr,1] <- b2$mean_SM05[which(b2$date == B2SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B2SIG0$xpol_fore[i] + B2SIG0$xpol_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B2SIG0$xpol_fore[i]
  SIG0vsSMC[cntr,3] <- B2SIG0$xpol_aft[i]
  SIG0vsSMC[cntr,4] <- "XPOL"
  SIG0vsSMC[cntr,5] <- "B1500"
  SIG0vsSMC[cntr,6] <- as.character.Date(b2$date[which(b2$date == B2SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B2SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B2SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B2SIG0$lia_aft[i]
  cntr <- cntr + 1

  
}

for (i in 1:length(B3SIG0$time)){
  
  #if (B3SIG0$dist[i] >= 0.006) {next}
  if (any(b3$date == B3SIG0$time[i]) == F) {next}
  #if (B3SIG0$rflag[i] >= 8 & B3SIG0$rflag[i] < 16) {next}
  #print(i)
  SIG0vsSMC[cntr,1] <- b3$mean_SM05[which(b3$date == B3SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B3SIG0$vv_fore[i] + B3SIG0$vv_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B3SIG0$vv_fore[i]
  SIG0vsSMC[cntr,3] <- B3SIG0$vv_aft[i]
  SIG0vsSMC[cntr,4] <- "VV"
  SIG0vsSMC[cntr,5] <- "B2000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b3$date[which(b3$date == B3SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B3SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B3SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B3SIG0$lia_aft[i]
  cntr <- cntr + 1

  SIG0vsSMC[cntr,1] <- b3$mean_SM05[which(b3$date == B3SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B3SIG0$hh_fore[i] + B3SIG0$hh_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B3SIG0$hh_fore[i]
  SIG0vsSMC[cntr,3] <- B3SIG0$hh_aft[i]
  SIG0vsSMC[cntr,4] <- "HH"
  SIG0vsSMC[cntr,5] <- "B2000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b3$date[which(b3$date == B3SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B3SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B3SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B3SIG0$lia_aft[i]
  cntr <- cntr + 1

  SIG0vsSMC[cntr,1] <- b3$mean_SM05[which(b3$date == B3SIG0$time[i])]
  #SIG0vsSMC[cntr,2] <- (B3SIG0$xpol_fore[i] + B3SIG0$xpol_aft[i]) / 2
  SIG0vsSMC[cntr,2] <- B3SIG0$xpol_fore[i]
  SIG0vsSMC[cntr,3] <- B3SIG0$xpol_aft[i]
  SIG0vsSMC[cntr,4] <- "XPOL"
  SIG0vsSMC[cntr,5] <- "B2000"
  SIG0vsSMC[cntr,6] <- as.character.Date(b3$date[which(b3$date == B3SIG0$time[i])])
  SIG0vsSMC[cntr,7] <- B3SIG0$nadir_dist[i]
  SIG0vsSMC[cntr,8] <- B3SIG0$lia_fore[i]
  SIG0vsSMC[cntr,9] <- B3SIG0$lia_aft[i]
  cntr <- cntr + 1

  
}

val <- which(SIG0vsSMC$sig0f != -9999 & SIG0vsSMC$sig0a != -9999)
SIG0vsSMC <- SIG0vsSMC[val,]


plotSIG0vsSMC <- function(lookdir, col, SIG0vsSMC){
  
  if (lookdir=="fwd"){
    
    rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
    r <- rf(32)
    
    #val <- which(SIG0vsSMC$sig0 != -9999)
    p1 <- ggplot(data=SIG0vsSMC, aes(x=smc,y=10*log10(sig0f) , colour=LIAf)) + 
          geom_point()  +
          theme(aspect.ratio=1) +
          geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
          #scale_colour_gradientn(colours=r) +
          facet_grid(. ~ pol) + ggtitle("in-situ SM vs SIG0 (Fwd)\n")
    
    dev.new()
    
    
    valb1 <- which(SIG0vsSMC$pol == "VV")
    valb2 <- which(SIG0vsSMC$pol == "HH")
    valb3 <- which(SIG0vsSMC$pol == "XPOL")
    
    print("VV")
    print(cor(SIG0vsSMC$smc[valb1], SIG0vsSMC$sig0f[valb1])^2)
    print("HH")
    print(cor(SIG0vsSMC$smc[valb2], SIG0vsSMC$sig0f[valb2])^2)
    print("XPOL")
    print(cor(SIG0vsSMC$smc[valb3], SIG0vsSMC$sig0f[valb3])^2)
    
    return(p1)
    
  }
  
  if (lookdir=="aft"){
    
    rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
    r <- rf(32)
    
    #val <- which(SIG0vsSMC$sig0 != -9999)
    p1 <- ggplot(data=SIG0vsSMC, aes(x=smc,y=10*log10(sig0a) , colour=LIAa)) + 
      geom_point()  +
      theme(aspect.ratio=1) +
      geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
      #scale_colour_gradientn(colours=r) +
      facet_grid(. ~ pol) + ggtitle("in-situ SM vs SIG0 (Aft)\n")
    
    dev.new()
    
    
    valb1 <- which(SIG0vsSMC$pol == "VV")
    valb2 <- which(SIG0vsSMC$pol == "HH")
    valb3 <- which(SIG0vsSMC$pol == "XPOL")
    
    print("VV")
    print(cor(SIG0vsSMC$smc[valb1], SIG0vsSMC$sig0a[valb1])^2)
    print("HH")
    print(cor(SIG0vsSMC$smc[valb2], SIG0vsSMC$sig0a[valb2])^2)
    print("XPOL")
    print(cor(SIG0vsSMC$smc[valb3], SIG0vsSMC$sig0a[valb3])^2)
    
    return(p1)
    
  }
  
  
}

plotNDISTvsLIA <- function(SIG0vsSMC){
  
  NDISTvsLIA <- data.frame(ndist=SIG0vsSMC$NDIST[SIG0vsSMC$pol == "VV"],
                           lia=SIG0vsSMC$LIAf[SIG0vsSMC$pol == "VV"],
                           stat=SIG0vsSMC$stat[SIG0vsSMC$pol == "VV"],
                           dir=rep("FWD", length(which(SIG0vsSMC$pol == "VV"))))
  NDISTvsLIA <- rbind(NDISTvsLIA,
                      data.frame(ndist=SIG0vsSMC$NDIST[SIG0vsSMC$pol == "VV"],
                                 lia=SIG0vsSMC$LIAa[SIG0vsSMC$pol == "VV"],
                                 stat=SIG0vsSMC$stat[SIG0vsSMC$pol == "VV"],
                                 dir=rep("AFT", length(which(SIG0vsSMC$pol == "VV")))))
  
  p1 <- ggplot(data=NDISTvsLIA, aes(x=ndist/1000,y=lia)) + 
    geom_point()  +
    theme(aspect.ratio=1) +
    geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
    facet_grid(stat ~ dir) + ggtitle("Nadir-Distance vs LIA\n")
  
  dev.new()
  return(p1)
  
}

plotLIAvsSIG0 <- function(lookdir, SIG0vsSMC){
  
  if (lookdir == "fwd"){
    p1 <- ggplot(data=SIG0vsSMC, aes(x=LIAf,y=10*log10(sig0f))) + 
      geom_point()  +
      theme(aspect.ratio=1) +
      geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
      facet_grid(stat ~ pol) + ggtitle("LIA vs SIG0 (Fwd)\n")
    
    dev.new()
  }
  
  if (lookdir == "aft"){
    p1 <- ggplot(data=SIG0vsSMC, aes(x=LIAa,y=10*log10(sig0a))) + 
      geom_point()  +
      theme(aspect.ratio=1) +
      geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
      facet_grid(stat ~ pol) + ggtitle("LIA vs SIG0 (Aft)\n")
    
    dev.new()
  }
  
  return(p1)
  
}




# 
# #SVR estimation test
# #reshaping sig0vsstation data frame
# SIG0vsSMC <- SIG0vsSMC[val,]
# 
# SVRdt <- data.frame(smc=SIG0vsSMC$smc[SIG0vsSMC$pol == "VV"],
#                     vv=SIG0vsSMC$sig0[SIG0vsSMC$pol == "VV"],
#                     hh=SIG0vsSMC$sig0[SIG0vsSMC$pol == "HH"],
#                     hv=SIG0vsSMC$sig0[SIG0vsSMC$pol == "XPOL"])
# 
# 
# #Splitting the dataset into training and test-set
# 
# nrofrows <- nrow(SVRdt)
# 
# #setting the size of the training set
# nTrain <- 23
# pTrain <- 23*(1/nrofrows)
# 
# #using a pseudo random sampling strategy
# set.seed(3456)
# trainIndices <- createDataPartition(SVRdt$smc, p=pTrain, list=FALSE, times=1)
# #load("./SMAP/trainingindices.dat")
# training <- SVRdt[trainIndices,]
# testing <- SVRdt[-trainIndices,]
# 
# #--------------------------------------------------------------------------------------------------
# #Building the SVR model
# #package e1071
# 
# #setting the tuning options
# SVRtuningsettings <- tune.control(sampling = "cross", 
#                                   sampling.aggregate = mean,
#                                   sampling.dispersion = sd, 
#                                   cross = 5,
#                                   best.model = TRUE,
#                                   performances = TRUE)
# 
# #tuning the model (finding the best set  hyper parameters)
# SVRtuning <- tune(svm, smc ~ vv + hh + hv, 
#                   data=training, 
#                   kernel="radial",
#                   ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
#                   tunecontrol=SVRtuningsettings)
# 
# tunedModel <- SVRtuning$best.model
# 
# #--------------------------------------------------------------------------------------------------
# #Estimating SMC
# 
# print("Overall performance based on full testset:")
# SMCpredicted <- predict(tunedModel, testing)
# error <- sqrt(mean((testing$smc - SMCpredicted)^2))
# r2 <- cor(testing$smc, y=SMCpredicted, method="pearson")^2
# print(paste("Error:",error,"R2:",r2))
# 
# tmp <- data.frame(true=testing$smc, predicted=SMCpredicted)
# 
# p2 <- ggplot(data=tmp, aes(x=true,y=predicted)) + 
#   geom_point()  +
#   theme(aspect.ratio=1) +
#   geom_abline(intercept = 0, slope = 1, linetype=2) + 
#   xlim(0,50) + ylim(0,50)
# 
# dev.new()
# p2
