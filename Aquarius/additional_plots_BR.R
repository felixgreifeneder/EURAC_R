#loading required packages
library("ggplot2")
library("reshape2")
library("hexbin")
library("RColorBrewer")
library("gridExtra")
library("raster")
library("R.matlab")

#load testset
load("X:/Workspaces/GrF/Processing/SCA_paper/era_valdiation_svr/Aquarius_with_means/train_test_mean_slope.dat")


#estimate SMC
matdata <- readMat("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/data&code for Bayes/file finali ASCAT/estim_smc_vv_f_vhvv_ts_n2.mat")
SMCpredicted <- drop(matdata[[1]])

#retrieve land-cover
globcover <- raster("X:/Workspaces/GrF/01_Data/ANCILLARY/GLOBCOVER/GLOBCOVER_L4_200901_200912_V2.3.tif")

t <- matrix(c(testset$lon, testset$lat), nrow(testset), 2)
lclist <- extract(globcover, t)
rslist <- data.frame(true=testset$smc, est=SMCpredicted, lc=as.factor(lclist), lat=testset$lat, lon=testset$lon)
uselc <- c(14, 20, 30, 40, 60, 90, 100, 110, 120, 130, 140, 150)

#rslist <- rslist[lclist %in% uselc,]

#create correlation per land-cover output ----------------------------------------------------------
outdir <- "X:/Workspaces/GrF/Processing/SCA_paper/ASCATvsAQU/"



# world map accuracy plot

# calculate rmse by location
uniquePos <- unique(rslist[,4:5])
nUL <- nrow(uniquePos)

errors <- data.frame(lat=rep(0, nUL), lon=rep(0, nUL), rmse=rep(0, nUL), r=rep(0, nUL))

for (i in c(1:nUL)){
  tlon <- uniquePos$lon[i]
  tlat <- uniquePos$lat[i]
  
  tr <- rslist$true[which((rslist$lon == tlon) & (rslist$lat == tlat))]
  est <- rslist$est[which((rslist$lon == tlon) & (rslist$lat == tlat))]
  
  error <- RMSE(est, tr)
  r <- cor(est, tr)
  errors$lat[i] <- tlat
  errors$lon[i] <- tlon
  errors$rmse[i] <- error
  errors$r[i] <- r
}

map <- borders("world", colour="gray50", fill="gray50")
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

path <- paste(outdir, 'br_rmse.png', sep="")

#dev.new(width=10, height=6, noRStudioGD=T)
png(filename = path, res=300, width = 2000, height = 1500)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=rmse), data=errors, size=0.8) + 
  scale_colour_gradientn(colours=r, name="RMSE", limit=c(0,0.3), na.value="red", breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("BR prediction - RMSE\n")

dev.off()

path <- paste(outdir, 'br_r.png', sep="")

#dev.new(width=10, height=6, noRStudioGD=T)
png(filename = path, res=300, width = 2000, height = 1500)
map <- borders("world", colour="gray50", fill="gray50")
#dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=r), data=errors, size=0.8) + 
  scale_colour_gradientn(colours=rev(r), name="R", labels=c("-1.0 ","-0.5","0.0","0.5","1.0"), limit=c(-1,1))
gridPlt + ggtitle("BR prediction - R\n")

dev.off()
