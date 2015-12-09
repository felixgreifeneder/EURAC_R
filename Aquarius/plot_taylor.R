library(plotrix)
library(RColorBrewer)

load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat/tss.dat")
ascat <- allData
load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_slope/tss.dat")
ascat_slope <- allData
load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_aqu/tss.dat")
ascat_aqu <- allData

ascat$lc <- as.factor(ascat$lc)
ascat_slope$lc <- as.factor(ascat_slope$lc)
ascat_aqu$lc <- as.factor(ascat_aqu$lc)

modelList <- list(ascat, ascat_slope, ascat_aqu)
levelNames <- levels(ascat$lc)
levelNames <- levelNames[1:11]

lcColours <- brewer.pal(11, "Paired")
dev.new(width=6, height=6, noRStudioGD=T)
taylor.diagram(ascat$smcts.insitu.mean, ascat$smcts.pred, pos.cor = F, col="white")

for (i in 1:3){
  for (lcN in 1:length(levelNames)){
    taylor.diagram(modelList[[i]]$smcts.insitu.mean[which(modelList[[i]]$lc == levelNames[lcN])], 
                   modelList[[i]]$smcts.pred[which(modelList[[i]]$lc == levelNames[lcN])], 
                   add=T, 
                   col=lcColours[lcN],
                   pch=i+14,
                   pos.cor = F,
                   pcex=1)
  }
}


