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

#load testset
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius/parameters_global_R")
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/train_test_mean_slope.dat")

#load SVR model
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/ascat_aqu/SVRmodel_performance.dat")

#estimate SMC
tunedModel <- SVRtuning$best.model
SMCpredicted <- parallel_predictions(tunedModel, testset)

#retrieve land-cover
globcover <- raster("X:/Workspaces/GrF/01_Data/ANCILLARY/GLOBCOVER/GLOBCOVER_L4_200901_200912_V2.3.tif")

t <- matrix(c(testset$lon, testset$lat), nrow(testset), 2)
lclist <- extract(globcover, t)
rslist <- data.frame(true=testset$smc, est=SMCpredicted, lc=as.factor(lclist))
uselc <- c(14, 20, 30, 40, 60, 90, 100, 110, 120, 130, 140, 150)

rslist <- rslist[lclist %in% uselc,]

#create correlation per land-cover output ----------------------------------------------------------
outdir <- "X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/ascat_aqu/"


sink(paste(outdir, "lc_metrics.txt", sep=""))

cat("Number of Samples: ")
cat(nrow(rslist))
cat("\n\n")
cat("Number of samples by LC \n")
for (i in uselc){
  cat(i)
  cat(": ")
  cat(length(which(rslist$lc == i)))
  cat("\n")
}
cat("\n")

cat("Overall correlation")
cat("\n")
cat(cor(rslist$true, rslist$est))
cat("\n")
cat("Correlations by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  cat(cor(rslist$true[userows], rslist$est[userows]))
  cat("\n")
}
cat("\n")
cat("Overall mean bias")
cat("\n")
cat(mean(rslist$est-rslist$true, na.rm=T))
cat("\n")
cat("Mean Bias by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  cat(mean(rslist$est[userows]-rslist$true[userows], na.rm=T))
  cat("\n")
}
cat("\n")
cat("Overall RMSE")
cat("\n")
cat(sqrt(mean((rslist$est-rslist$true)^2, na.rm=T)))
cat("\n")
cat("RMSE by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  cat(sqrt(mean((rslist$est[userows]-rslist$true[userows])^2, na.rm=T)))
  cat("\n")
}
cat("\n")

cat("Overall MSE")
cat("\n")
cat(mean((rslist$est-rslist$true)^2, na.rm=T))
cat("\n")
cat("RMSE by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  cat(mean((rslist$est[userows]-rslist$true[userows])^2, na.rm=T))
  cat("\n")
}
cat("\n")

cat("Overall Slope")
cat("\n")
tvs.lm <- lm(true ~ est, data=rslist)
cat(coefficients(tvs.lm)[2])
cat("\n")
cat("Slope by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  tvs.lm <- lm(true ~ est, data=rslist[userows,])
  cat(coefficients(tvs.lm)[2])
  cat("\n")
}
cat("\n")
cat("Overall Intercept")
cat("\n")
tvs.lm <- lm(true ~ est, data=rslist)
cat(coefficients(tvs.lm)[1])
cat("\n")
cat("Intercept by LC")
cat("\n")
for (i in uselc){
  userows <- which(rslist$lc == i)
  cat(i)
  cat(": ")
  tvs.lm <- lm(true ~ est, data=rslist[userows,])
  cat(coefficients(tvs.lm)[1])
  cat("\n")
}
cat("\n")



sink()

#plot ------------------------------------------------
#histo scatterplot
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)


#for (i in levels(rslist$lc)){
  
  #path <- paste(outdir, 'denspl/', i, '.png', sep="")
path <- paste(outdir, 'truevsest_dens.png', sep="")

  #vlc <- which(lclist == 14 | lclist == 20 | lclist == 30 | lclist == 60 | lclist ==110 | lclist == 120 |
  #               lclist == 130 | lclist == 140 | lclist == 150)
  #vlc <- which(lclist == 30)
  #vlc <- which(lclist == i)

  rslist_flt <- rslist#[vlc,]

  set.seed(3456)
  #testIndices <- createDataPartition(rslist_flt$true, p=0.01, list=FALSE, times=1)
  #rslist_flt <- rslist_flt[testIndices,]

  #dev.new(width=6, height=6, noRStudioGD = T)
  png(filename = path, res=300, width = 1500, height = 1500)
  
  empty <- ggplot() + geom_point(aes(1,1), colour="white") +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank()
    )

  scatter <- ggplot(rslist_flt, aes(true, est)) + stat_bin2d(bins=100) + geom_smooth(method=lm, se=F, fullrange=T) +
    xlab("True SMC [m3m-3]") + ylab("Estimated SMC [m3m-3]") + xlim(0,0.5) + ylim(0,0.5) +
    theme(aspect.ratio=1, legend.position="none") + scale_fill_gradientn(colours=r, trans="log")

  plot_top <- ggplot(rslist_flt, aes(true)) + 
    geom_density(alpha=.5) + theme(axis.title.x = element_blank()) + xlim(0,0.5)
    theme(legend.position = "none")

  plot_right <- ggplot(rslist_flt, aes(est)) + 
    geom_density(alpha=.5) + 
    coord_flip() + theme(axis.title.y = element_blank()) + xlim(0,0.5)
    theme(legend.position = "none") 

  grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  #ggsave(filename=path)
  dev.off()
  
#}
  
  

