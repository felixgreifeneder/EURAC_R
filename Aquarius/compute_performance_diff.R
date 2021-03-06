library("psych")
library("scales")
library("ggplot2")
library("RColorBrewer")
library("gridExtra")
library("R.matlab")
library("raster")
source("./Aquarius/parallel_predictions.R")

#compute difference between different estimation approaches
#additionally compute significance measures

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/train_test_mean_slope.dat")







#get performance measures ASCAT
#--------------------------------------------------------------------------------------------------

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/ascat/SVRmodel_performance.dat")
tunedModel <- SVRtuning$best.model

print("Overall performance based on full testset - ASCAT:")
estimates <- parallel_predictions(tunedModel, testset)
error <- sqrt(mean((testset$smc - estimates)^2))
r2 <- cor(testset$smc, y=estimates, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

ascat <- prediction_performance


#get performance measures ASCAT+slope
#--------------------------------------------------------------------------------------------------

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/ascat_slope/SVRmodel_performance.dat")
tunedModel <- SVRtuning$best.model

print("Overall performance based on full testset - ASCAT+slope:")
estimates <- parallel_predictions(tunedModel, testset)
error <- sqrt(mean((testset$smc - estimates)^2))
r2 <- cor(testset$smc, y=estimates, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

ascat_slope <- prediction_performance


#get performance measures ASCAT
#--------------------------------------------------------------------------------------------------

load("C:/Users/FGreifeneder/Documents/tmp_proc/Aquarius_with_means/ascat_aqu/SVRmodel_performance.dat")
tunedModel <- SVRtuning$best.model

print("Overall performance based on full testset - ASCAT+aqu:")
estimates <- parallel_predictions(tunedModel, testset)
error <- sqrt(mean((testset$smc - estimates)^2))
r2 <- cor(testset$smc, y=estimates, method="pearson")^2
print(paste("Error:",error,"R2:",r2))

ascat_aqu <- prediction_performance


#--------------------------------------------------------------------------------------------------------
#compute diff 
#--------------------------------------------------------------------------------------------------------

globcover <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/GLOBCOVER_L4_200901_200912_V2.3.tif")

nlocations <- nrow(ascat)
diff_data <- data.frame(lat=rep(-9999,nlocations), 
                        lon=rep(-9999,nlocations),
                        rdiff_slope=rep(-9999,nlocations),
                        rdiff_aqu=rep(-9999,nlocations),
                        rmsediff_slope=rep(-9999,nlocations),
                        rmsediff_aqu=rep(-9999,nlocations),
                        msediff_slope=rep(-9999,nlocations),
                        msediff_aqu=rep(-9999,nlocations),
                        rconf_l=rep(-9999,nlocations),
                        rconf_u=rep(-9999,nlocations),
                        rconf_span = rep(-9999,nlocations),
                        rftest_slope = rep(-9999,nlocations),
                        rftest_aqu = rep(-9999,nlocations),
                        tt_slope = rep(-9999,nlocations),
                        tt_aqu = rep(-9999,nlocations),
                        tv_slope = rep(-9999,nlocations),
                        tv_aqu = rep(-9999,nlocations),
                        lc = rep(-9999,nlocations))

#for each location compute different metrics
for (i in c(1:nlocations)){
  if (ascat$n[i] <= 3) next
  
  diff_data$lat[i] <- ascat$lat[i]
  diff_data$lon[i] <- ascat$lon[i]
  
  diff_data$rdiff_slope[i] <- ascat_slope$r2[i] - ascat$r2[i]
  diff_data$rdiff_aqu[i] <- ascat_aqu$r2[i] - ascat$r2[i]
      
  rconf <- r.con(sqrt(ascat$r2[i]), ascat$n[i], p=.95)
  rftest_slope <- (fisherz(sqrt(ascat$r2[i]))-fisherz(sqrt(ascat_slope$r2[i])))/sqrt(1/(ascat$n[i]-3)+1/(ascat_slope$n[i]-3))
  rftest_aqu <- (fisherz(sqrt(ascat$r2[i]))-fisherz(sqrt(ascat_aqu$r2[i])))/sqrt(1/(ascat$n[i]-3)+1/(ascat_aqu$n[i]-3))
  diff_data$rconf_l[i] <- rconf[1]
  diff_data$rconf_u[i] <- rconf[2]
  diff_data$rconf_span[i] <- rconf[2] - rconf[1]
  diff_data$rftest_slope[i] <- rftest_slope
  diff_data$rftest_aqu[i] <- rftest_aqu
  
  diff_data$rmsediff_slope[i] <- ascat_slope$rmse[i] - ascat$rmse[i]
  diff_data$rmsediff_aqu[i] <- ascat_aqu$rmse[i] - ascat$rmse[i]
  tt <- (ascat_slope$mse[i]-ascat$mse[i])/sqrt(ascat_slope$vari[i]/ascat_slope$n[i] + ascat$vari[i]/ascat$n[i])
  v <- (ascat_slope$vari[i]/ascat_slope$n[i] + ascat$vari[i]/ascat$n[i])^2/
    ((ascat_slope$vari[i]/ascat_slope$n[i])^2/(ascat_slope$n[i]-1) + (ascat$vari[i]/ascat$n[i])^2/(ascat$n[i]-1))
  diff_data$tt_slope[i] <- tt
  diff_data$tv_slope[i] <- qt(1-0.025, df=round(v))
  tt <- (ascat_aqu$mse[i]-ascat$mse[i])/sqrt(ascat_aqu$vari[i]/ascat_aqu$n[i] + ascat$vari[i]/ascat$n[i])
  v <- (ascat_aqu$vari[i]/ascat_aqu$n[i] + ascat$vari[i]/ascat$n[i])^2/
    ((ascat_aqu$vari[i]/ascat_aqu$n[i])^2/(ascat_aqu$n[i]-1) + (ascat$vari[i]/ascat$n[i])^2/(ascat$n[i]-1))
  diff_data$tt_aqu[i] <- tt
  diff_data$tv_aqu[i] <- qt(1-0.025, df=round(v))
  diff_data$lc[i] <- as.numeric(extract(globcover, matrix(c(ascat$lon[i], ascat$lat[i]), 1, 2)))
  
}

val <- which(diff_data$lat != -9999)

#normalize differences
diff_data$rdiff_slope[val] <- diff_data$rdiff_slope[val] / abs(max(diff_data$rdiff_slope[val], na.rm=T))
diff_data$rdiff_aqu[val] <- diff_data$rdiff_aqu[val] / abs(max(diff_data$rdiff_aqu[val], na.rm=T))
diff_data$rmsediff_slope[val] <- diff_data$rmsediff_slope[val] / abs(max(diff_data$rmsediff_slope[val], na.rm=T))
diff_data$rmsediff_aqu[val] <- diff_data$rmsediff_aqu[val] / abs(max(diff_data$rmsediff_aqu[val], na.rm=T))

rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

#correlation change significance ascat - ascat+slope
#---------------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, z=abs(diff_data$rftest_slope), issignif=as.factor(abs(diff_data$rftest_slope) > 1.96))

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=issignif), data=tmp[val,], size=1.2) + 
  scale_colour_discrete(name="z")
gridPlt + ggtitle("Significant change R: ASCAT - ASCAT+slope\n")

#correlation change significance ascat - ascat+aqu
#-------------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, z=abs(diff_data$rftest_aqu), issignif=as.factor(abs(diff_data$rftest_aqu) > 1.96))

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=issignif), data=tmp[val,], size=1.2) + 
  scale_colour_discrete(name="z")
gridPlt + ggtitle("Significant change R: ASCAT - ASCAT+VH/VV\n")

#rmse change significance ascat - ascat+slope
#--------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, t=abs(diff_data$tt_slope), issignif=as.factor(abs(diff_data$tt_slope) > diff_data$tv_slope))

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=issignif), data=tmp[val,], size=1.2) + 
  scale_colour_discrete(name="z")
gridPlt + ggtitle("Significant change RMSE: ASCAT - ASCAT+slope\n")

#rmse change significance ascat - ascat+aqu
#--------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, t=abs(diff_data$tt_aqu), issignif=as.factor(abs(diff_data$tt_aqu) > diff_data$tv_aqu))

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=issignif), data=tmp[val,], size=1.2) + 
  scale_colour_discrete(name="z")
gridPlt + ggtitle("Significant change RMSE: ASCAT - ASCAT+aqu\n")

#rmse diff
#---------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, diff=diff_data$rmsediff_slope)
tmp$diff[which(tmp$diff > 0.25)] <- 0.25
tmp$diff[which(tmp$diff < -0.25)] <- -0.25
  
map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=diff), data=tmp[val,], size=1.2) + 
  scale_colour_gradientn(colours=r, name="dRMSE", limit=c(-0.25,0.25), 
                         breaks=c(-0.25,-0.125,0,0.125,0.25), 
                         labels=c("< -0.25", "-0.125", "0", "0.125", "> 0.25"))
gridPlt + ggtitle("RMSE difference: ASCAT+slope - ASCAT\n")

tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, diff=diff_data$rmsediff_aqu)
tmp$diff[which(tmp$diff > 0.25)] <- 0.25
tmp$diff[which(tmp$diff < -0.25)] <- -0.25

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=diff), data=tmp[val,], size=1.2) + 
  scale_colour_gradientn(colours=r, name="dRMSE", limit=c(-0.25,0.25), 
                         breaks=c(-0.25,-0.125,0,0.125,0.25), 
                         labels=c("< -0.25", "-0.125", "0", "0.125", "> 0.25"))
gridPlt + ggtitle("RMSE difference: ASCAT+VH/VV - ASCAT\n")

#rdiff
----------------------------------------------
tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, diff=diff_data$rdiff_slope)
tmp$diff[which(tmp$diff > 0.25)] <- 0.25
tmp$diff[which(tmp$diff < -0.25)] <- -0.25

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=diff), data=tmp[val,], size=1.2) + 
  scale_colour_gradientn(colours=rev(r), name="dR", limit=c(-0.25,0.25), 
                         breaks=c(-0.25,-0.125,0,0.125,0.25), 
                         labels=c("< -0.25", "-0.125", "0", "0.125", "> 0.25"))
gridPlt + ggtitle("R2 difference: ASCAT+slope - ASCAT\n")

tmp <- data.frame(lat=diff_data$lat, lon=diff_data$lon, diff=diff_data$rdiff_aqu)
tmp$diff[which(tmp$diff > 0.25)] <- 0.25
tmp$diff[which(tmp$diff < -0.25)] <- -0.25

map <- borders("world", colour="gray50", fill="gray50")
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=diff), data=tmp[val,], size=1.2) + 
  scale_colour_gradientn(colours=rev(r), name="dR", limit=c(-0.25,0.25), 
                         breaks=c(-0.25,-0.125,0,0.125,0.25), 
                         labels=c("< -0.25", "-0.125", "0", "0.125", "> 0.25"))
gridPlt + ggtitle("R2 difference: ASCAT+VH/VV - ASCAT\n")

#boxplots
#-----------------------------------------------------------------------------

diff_data$lc <- as.factor(diff_data$lc)
new_lc <- diff_data$lc
new_lc <- as.character(new_lc)
new_lc[which(diff_data$lc == 14 | diff_data$lc == 20 | diff_data$lc == 30)] <- "croplands"
new_lc[which(diff_data$lc == 40 | diff_data$lc == 60 | diff_data$lc == 90 | diff_data$lc == 100)] <- "forest, open"
new_lc[which(diff_data$lc == 50 | diff_data$lc == 70)] <- "forest, closed"
new_lc[which(diff_data$lc == 110 | diff_data$lc == 120 | diff_data$lc == 130 | diff_data$lc == 140)] <- "shrublands/grasslands"
new_lc[which(diff_data$lc == 150 | diff_data$lc == 200)] <- "sparsely vegetated/bare"
new_lc[which(diff_data$lc == 160 | diff_data$lc == 170 | diff_data$lc == 180)] <- "regularly flooded"
new_lc[which(diff_data$lc == 190)] <- "urban"

diff_data <- data.frame(diff_data, lc_labels=new_lc)

use_lc <- which((diff_data$lc == 14 | diff_data$lc == 20 | diff_data$lc == 30 | diff_data$lc == 40 | 
                   diff_data$lc == 60 | diff_data$lc == 90 | diff_data$lc == 100 | diff_data$lc == 50 | 
                   diff_data$lc == 70 | diff_data$lc == 110 | diff_data$lc == 120 | diff_data$lc == 130 | 
                   diff_data$lc == 140 | diff_data$lc == 150 | diff_data$lc == 160 | diff_data$lc == 170 | 
                   diff_data$lc == 200 | diff_data$lc == 190 | diff_data$lc == 180) & 
                  diff_data$lat != -9999)
val <- use_lc

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(diff_data[val,], aes(x=lc_labels, y=rdiff_slope, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(-0.5,0.5) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("R2 difference") + 
  ggtitle("R2 difference: ASCAT+slope - ASCAT\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(diff_data[val,], aes(x=lc_labels, y=rdiff_aqu, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(-0.5,0.5) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("R2 difference") + 
  ggtitle("R2 difference: ASCAT+VH/VV - ASCAT\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(diff_data[val,], aes(x=lc_labels, y=rmsediff_slope, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(-0.5,0.5) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("RMSE difference") + 
  ggtitle("RMSE difference: ASCAT+slope - ASCAT\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))

dev.new(width=8, height=6, noRStudioGD=T)
ggplot(diff_data[val,], aes(x=lc_labels, y=rmsediff_aqu, fill=lc)) + geom_boxplot() + 
  guides(fill=FALSE) + 
  ylim(-0.5,0.5) + 
  xlab("GLOBCOVER land-cover class") + 
  ylab("RMSE difference") + 
  ggtitle("RMSE difference: ASCAT+VH/VV - ASCAT\n") +
  theme(axis.text.x  = element_text(angle=45, vjust=0.5))