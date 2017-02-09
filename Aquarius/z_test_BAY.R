#loading required packages
library("e1071")
library("caret")
library("ggplot2")
library("reshape2")
library("hexbin")
library("RColorBrewer")
library("gridExtra")
source("./Aquarius/parallel_predictions.R")
library("raster")
library("gridExtra")
library("psych")

load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/train_test_mean_slope.dat")

matdata <- readMat("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/data&code for Bayes/file finali ASCAT/estim_smc_vv_f_ts_n.mat")
SMCpredicted_ascat <- drop(matdata[[1]])
matdata <- readMat("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/data&code for Bayes/file finali ASCAT/estim_smc_vv_f_slope_ts_n2.mat")
SMCpredicted_slope <- drop(matdata[[1]])
matdata <- readMat("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/data&code for Bayes/file finali ASCAT/estim_smc_vv_f_vhvv_ts_n2.mat")
SMCpredicted_aqu <- drop(matdata[[1]])


#retrieve land-cover
globcover <- raster("X:/Workspaces/GrF/01_Data/ANCILLARY/GLOBCOVER/GLOBCOVER_L4_200901_200912_V2.3.tif")

t <- matrix(c(testset$lon, testset$lat), nrow(testset), 2)
lclist <- extract(globcover, t)
rslist <- data.frame(true=testset$smc, 
                     ascat=SMCpredicted_ascat, 
                     slope=SMCpredicted_slope,
                     aqu=SMCpredicted_aqu,
                     lc=as.factor(lclist))
uselc <- c(14, 20, 30, 40, 60, 90, 100, 110, 120, 130, 140, 150)

rslist <- rslist[lclist %in% uselc,]

#create correlation per land-cover output ----------------------------------------------------------
outdir <- "X:/Workspaces/GrF/Processing/SCA_paper/era_validation_bay/"

sink(paste(outdir, "lc_change.txt", sep=""))
cat("Z-Test for significance of correlation change, by land cover: \n\n")
cat("LC          SLOPE-ASCAT          AQU-ASCAT\n")
cat("------------------------------------------\n\n")

for (i in uselc){
  
  cat(i)
  cat("           ")
  
  #test difference between r ascat and r ascat slope
  
  r1 <- cor(rslist$true[which(rslist$lc==i)], rslist$ascat[which(rslist$lc==i)])
  r2 <- cor(rslist$true[which(rslist$lc==i)], rslist$slope[which(rslist$lc==i)])
  r3 <- cor(rslist$ascat[which(rslist$lc==i)], rslist$slope[which(rslist$lc==i)])
  r2_1 <- r1^2
  r2_2 <- r2^2
  r2_3 <- r3^2
  
  #n <- nrow(rslist)
  n <- length(which(rslist$lc==i))
  df <- n-3
  
  Zr1 <- fisherz(r1)
  Zr2 <- fisherz(r2)
  Zmr <- 0.5*(Zr1+Zr2)
  rm <- fisherz2r(Zmr)
  r2m <- rm^2
  
  s <- (r3*(1-r2m-r2m)-0.5*r2m*(1-r2m-r2m-r2_3))/((1-r2m)*(1-r2m))
  ZH <- (sqrt(df)*(Zr1-Zr2))/sqrt(2-2*s)
  
  p_sl <- 2*pnorm(-abs(ZH))
  
  cat(p_sl)
  cat("            ")
  
  #test difference between r ascat and r ascat aqu
  
  r1 <- cor(rslist$true[which(rslist$lc==i)], rslist$ascat[which(rslist$lc==i)])
  r2 <- cor(rslist$true[which(rslist$lc==i)], rslist$aqu[which(rslist$lc==i)])
  r3 <- cor(rslist$ascat[which(rslist$lc==i)], rslist$aqu[which(rslist$lc==i)])
  r2_1 <- r1^2
  r2_2 <- r2^2
  r2_3 <- r3^2
  
  n <- length(which(rslist$lc==i))
  df <- n-3
  
  Zr1 <- fisherz(r1)
  Zr2 <- fisherz(r2)
  Zmr <- 0.5*(Zr1+Zr2)
  rm <- fisherz2r(Zmr)
  r2m <- rm^2
  
  s <- (r3*(1-r2m-r2m)-0.5*r2m*(1-r2m-r2m-r2_3))/((1-r2m)*(1-r2m))
  ZH <- (sqrt(df)*(Zr1-Zr2))/sqrt(2-2*s)
  
  p_aq <- 2*pnorm(-abs(ZH))
  
  cat(p_aq)
  cat("\n")
  
}

cat("\n")

cat("T-Test for significance of Error improvement, by land cover: \n\n")
cat("LC          SLOPE-ASCAT          AQU-ASCAT\n")
cat("------------------------------------------\n\n")


for (i in uselc){
  
  cat(i)
  cat("           ")
  
  dev_ascat <- rslist$true[which(rslist$lc==i)] - rslist$ascat[which(rslist$lc==i)]
  dev_slope <- rslist$true[which(rslist$lc==i)] - rslist$slope[which(rslist$lc==i)]
  dev_aqu <- rslist$true[which(rslist$lc==i)] - rslist$aqu[which(rslist$lc==i)]
  
  tmp <- t.test(dev_ascat, dev_slope)
  tp_slp <- tmp$p.value
  
  tmp <- t.test(dev_ascat, dev_aqu)
  tp_aqu <- tmp$p.value
  
  cat(tp_slp)
  cat("            ")
  cat(tp_aqu)
  cat("\n")
  
  
  
}

sink()
