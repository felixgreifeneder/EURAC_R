#read data

outpath <- "C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_diffs_fltd/"

#ascat
load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_fltd/tss.dat")
#filter
diff = allData$smcts.era - allData$smcts.insitu.mean
v_era <- which(abs(diff) < 0.05 & allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210)
allData <- allData[v_era,]
#rename
ascat <- allData
ascat$lc <- as.factor(ascat$lc)
rm(allData)

#ascat+slope
load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_slope_fltd/tss.dat")
#filter
diff = allData$smcts.era - allData$smcts.insitu.mean
v_era <- which(abs(diff) < 0.05 & allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210)
allData <- allData[v_era,]
#rename
ascat_slope <- allData
ascat_slope$lc <- as.factor(ascat_slope$lc)
rm(allData)

#ascat+aqu
load("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_validation/svr_ascat_aqu_fltd/tss.dat")
#filter
diff = allData$smcts.era - allData$smcts.insitu.mean
v_era <- which(abs(diff) < 0.05 & allData$lc != 70 & allData$lc != 50 & allData$lc != 220 & allData$lc != 210)
allData <- allData[v_era,]
#rename
ascat_aqu <- allData
ascat_aqu$lc <- as.factor(ascat_aqu$lc)
rm(allData)


#compute performance metrics
#--------------------------------------------------------------------------------------------------
#bias
ascat_bias <- ascat$smcts.pred - ascat$smcts.insitu.mean
ascat_slope_bias <- ascat_slope$smcts.pred - ascat_slope$smcts.insitu.mean
ascat_aqu_bias <- ascat_aqu$smcts.pred - ascat_aqu$smcts.insitu.mean

#rmse and correlation by landcover
ascat_rmse <- rep(0, length(ascat_bias))
ascat_slope_rmse <- rep(0, length(ascat_bias))
ascat_aqu_rmse <- rep(0, length(ascat_bias))
ascat_r <- rep(0, length(ascat_bias))
ascat_slope_r <- rep(0, length(ascat_bias))
ascat_aqu_r <- rep(0, length(ascat_bias))


for (i in levels(ascat$lc)){
  ids <- which(ascat$lc == i)
  ascat_rmse[ids] <- sqrt(mean(ascat_bias[ids]^2, na.rm=T))
  ascat_slope_rmse[ids] <- sqrt(mean(ascat_slope_bias[ids]^2, na.rm=T))
  ascat_aqu_rmse[ids] <- sqrt(mean(ascat_aqu_bias[ids]^2, na.rm=T))
  ascat_r[ids] <- cor(ascat$smcts.pred[ids], ascat$smcts.insitu.mean[ids], use="pairwise.complete.obs")
  ascat_slope_r[ids] <- cor(ascat_slope$smcts.pred[ids], ascat_slope$smcts.insitu.mean[ids], use="pairwise.complete.obs")
  ascat_aqu_r[ids] <- cor(ascat_aqu$smcts.pred[ids], ascat_aqu$smcts.insitu.mean[ids], use="pairwise.complete.obs")
}

#combine dataframes
ascat$bias <- ascat_bias
ascat <- data.frame(ascat,
                    rmse = ascat_rmse,
                    r = ascat_r)
ascat_slope$bias <- ascat_slope_bias
ascat_slope <- data.frame(ascat_slope,
                          rmse = ascat_slope_rmse,
                          r = ascat_slope_r)
ascat_aqu$bias <- ascat_aqu_bias
ascat_aqu <- data.frame(ascat_aqu,
                        rmse = ascat_aqu_rmse,
                        r = ascat_aqu_r)

#compute improvements
#--------------------------------------------------------------------------------------------------

diff_ascat_slope <- data.frame(dbias = ascat_slope$bias - ascat$bias,
                               drmse = ascat_slope$rmse - ascat$rmse,
                               dr = ascat_slope$r - ascat$r)
#normalize
diff_ascat_slope <- data.frame(dbias = diff_ascat_slope$dbias/max(abs(diff_ascat_slope$dbias),na.rm=T),
                               drmse = diff_ascat_slope$drmse/max(abs(diff_ascat_slope$drmse),na.rm=T),
                               dr = diff_ascat_slope$dr/max(abs(diff_ascat_slope$dr),na.rm=T),
                               lc = ascat$lc)

diff_ascat_aqu <- data.frame(dbias = ascat_aqu$bias - ascat$bias,
                               drmse = ascat_aqu$rmse - ascat$rmse,
                               dr = ascat_aqu$r - ascat$r)

#normalize
diff_ascat_aqu <- data.frame(dbias = diff_ascat_aqu$dbias/max(abs(diff_ascat_aqu$dbias),na.rm=T),
                               drmse = diff_ascat_aqu$drmse/max(abs(diff_ascat_aqu$drmse),na.rm=T),
                               dr = diff_ascat_aqu$dr/max(abs(diff_ascat_aqu$dr),na.rm=T),
                               lc = ascat$lc)

#plots
#--------------------------------------------------------------------------------------------------

#bias boxplots
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=dbias, fill=lc), data=diff_ascat_slope[which(diff_ascat_slope$lc != 220 & diff_ascat_slope$lc != 210),]) + 
  geom_boxplot() + xlab("\nLand-Cover") + ylab("normalized diff. bias\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_dbias_ascat_slope.png", sep=""))


dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=dbias, fill=lc), data=diff_ascat_aqu[which(diff_ascat_aqu$lc != 220 & diff_ascat_aqu$lc != 210),]) + 
  geom_boxplot() + xlab("\nLand-Cover") + ylab("normalized diff. bias\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_dbias_ascat_aqu.png", sep=""))

#rmse barplot
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=drmse, fill=lc), data=diff_ascat_slope[which(diff_ascat_slope$lc != 220 & diff_ascat_slope$lc != 210),]) + 
  geom_bar(stat="identity") + xlab("\nLand-Cover") + ylab("normalized diff. rmse\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_drmse_ascat_slope.png", sep=""))

dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=drmse, fill=lc), data=diff_ascat_aqu[which(diff_ascat_aqu$lc != 220 & diff_ascat_aqu$lc != 210),]) + 
  geom_bar(stat="identity") + xlab("\nLand-Cover") + ylab("normalized diff. rmse\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_drmse_ascat_aqu.png", sep=""))

#r barplot
dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=dr, fill=lc), data=diff_ascat_slope[which(diff_ascat_slope$lc != 220 & diff_ascat_slope$lc != 210),]) + 
  geom_bar(stat="identity") + xlab("\nLand-Cover") + ylab("normalized diff. R\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_dr_ascat_slope.png", sep=""))

dev.new(width=7, height=5, noRStudioGD=T)
ggplot(aes(x=lc, y=dr, fill=lc), data=diff_ascat_aqu[which(diff_ascat_aqu$lc != 220 & diff_ascat_aqu$lc != 210),]) + 
  geom_bar(stat="identity") + xlab("\nLand-Cover") + ylab("normalized diff. R\n") + ylim(-1,1) +
  scale_fill_discrete(name="Land-Cover\n", 
                      labels=c("Cropland (50-70%)/vegetation (20-50%)",
                               "Vegetation (50-70%)/cropland (50-70%)",
                               "Closed to open mixed forest",
                               "Closed to open shrubland",
                               "Closed to open grassland",
                               "Bare areas"))

ggsave(filename=paste(outpath, "boxplot_dr_ascat_aqu.png", sep=""))

#output metrics to textfile

sink(paste(outpath, "metrics.txt", sep=""))

cat("Correlation improvements between ASCAT and ASCAT+slope:")
cat("\n")
for (i in levels(diff_ascat_slope$lc)){
  
  ids <- which(diff_ascat_slope$lc == i)
  cat(i)
  cat(": ")
  cat(diff_ascat_slope$dr[ids[1]])
  cat("\n")
  
}

cat("\n")
cat("RMSE imporvements between ASCAT and ASCAT+slope:")
cat("\n")
for (i in levels(diff_ascat_slope$lc)){
  
  ids <- which(diff_ascat_slope$lc == i)
  cat(i)
  cat(": ")
  cat(diff_ascat_slope$drmse[ids[1]])
  cat("\n")
  
}
cat("\n")
cat("\n")
cat("Correlation improvements between ASCAT and ASCAT+aqu:")
cat("\n")
for (i in levels(diff_ascat_aqu$lc)){
  
  ids <- which(diff_ascat_aqu$lc == i)
  cat(i)
  cat(": ")
  cat(diff_ascat_aqu$dr[ids[1]])
  cat("\n")
  
}

cat("\n")
cat("RMSE imporvements between ASCAT and ASCAT+aqu:")
cat("\n")
for (i in levels(diff_ascat_aqu$lc)){
  
  ids <- which(diff_ascat_aqu$lc == i)
  cat(i)
  cat(": ")
  cat(diff_ascat_aqu$drmse[ids[1]])
  cat("\n")
  
}
cat("\n")
sink()


