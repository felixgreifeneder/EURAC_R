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
rslist <- data.frame(true=testset$smc, est=SMCpredicted, ascatM=testset$ascatM, aquVV=testset$aquVV,lc=as.factor(lclist), lon=testset$lon, lat=testset$lat)
uselc <- c(14, 20, 30, 40, 60, 90, 100, 110, 120, 130, 140, 150)

#rslist <- rslist[lclist %in% uselc,]

#create correlation per land-cover output ----------------------------------------------------------
outdir <- "X:/Workspaces/GrF/Processing/SCA_paper/ASCATvsAQU/"



#plot ------------------------------------------------
#histo scatterplot
rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

# ASCAT VS AQUARIUS

path <- paste(outdir, 'ascatfwd_vs_aquvv.png', sep="")

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

scatter <- ggplot(dataCollection, aes(ascatF, aquVV)) + stat_bin2d(bins=100) + geom_smooth(method=lm, se=F, fullrange=T) +
  xlab("ASCAT Fwd (VV) [dB]") + ylab("AQUARIUS Mid (VV) [dB]") + xlim(-30,0) + ylim(-30,0) +
  theme(aspect.ratio=1, legend.position="none") + scale_fill_gradientn(colours=r, trans="log")

plot_top <- ggplot(dataCollection, aes(ascatF)) + 
  geom_density(alpha=.5) + theme(axis.title.x = element_blank()) + xlim(-30,0) + 
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none")

plot_right <- ggplot(dataCollection, aes(aquVV)) + 
  geom_density(alpha=.5) + 
  coord_flip() + theme(axis.title.y = element_blank()) + xlim(-30,0) +
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none") 

grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
  #ggsave(filename=path)
dev.off()
  
# ASCAT FWD vs ASCAT Mid

path <- paste(outdir, 'ascatfwd_vs_ascatmid.png', sep="")

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

scatter <- ggplot(dataCollection, aes(ascatF, ascatM)) + stat_bin2d(bins=100) + geom_smooth(method=lm, se=F, fullrange=T) +
  xlab("ASCAT Fwd (VV) [dB]") + ylab("ASCAT Mid (VV) [dB]") + xlim(-30,0) + ylim(-30,0) +
  theme(aspect.ratio=1, legend.position="none") + scale_fill_gradientn(colours=r, trans="log")

plot_top <- ggplot(dataCollection, aes(ascatF)) + 
  geom_density(alpha=.5) + theme(axis.title.x = element_blank()) + xlim(-30,0) + 
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none")

plot_right <- ggplot(dataCollection, aes(ascatM)) + 
  geom_density(alpha=.5) + 
  coord_flip() + theme(axis.title.y = element_blank()) + xlim(-30,0) +
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none") 

grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
#ggsave(filename=path)
dev.off()

# SMC vs ASCAT

path <- paste(outdir, 'smc_vs_ascat.png', sep="")

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

scatter <- ggplot(dataCollection, aes(smc, ascatM)) + stat_bin2d(bins=100) + geom_smooth(method=gam, formula=y~log(x), se=F, fullrange=T) +
  xlab("ERA-Land SMC [m3m-3]") + ylab("ASCAT (VV) [dB]") + ylim(-25,0) +
  theme(aspect.ratio=1, legend.position="none") + scale_fill_gradientn(colours=r, trans="log")

plot_top <- ggplot(dataCollection, aes(smc)) + 
  geom_density(alpha=.5) + theme(axis.title.x = element_blank()) + #xlim(-30,0) + 
  scale_y_continuous(breaks=c(0,4), labels=c("0","4"), limits=c(0,4)) +
  theme(legend.position = "none")

plot_right <- ggplot(dataCollection, aes(ascatM)) + 
  geom_density(alpha=.5) + 
  coord_flip() + theme(axis.title.y = element_blank()) + #xlim(-30,0) +
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none") 

grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
#ggsave(filename=path)
dev.off()

# SMC vs AQUARIUS

path <- paste(outdir, 'smc_vs_aqu.png', sep="")

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

scatter <- ggplot(dataCollection, aes(smc, aquVV)) + stat_bin2d(bins=100) + geom_smooth(method=gam, formula=y~log(x), se=F, fullrange=T) +
  xlab("ERA-Land SMC [m3m-3]") + ylab("AQUARIUS (VV) [dB]") + ylim(-25,0) +
  theme(aspect.ratio=1, legend.position="none") + scale_fill_gradientn(colours=r, trans="log")

plot_top <- ggplot(dataCollection, aes(smc)) + 
  geom_density(alpha=.5) + theme(axis.title.x = element_blank()) + #xlim(-30,0) + 
  scale_y_continuous(breaks=c(0,4), labels=c("0","4"), limits=c(0,4)) +
  theme(legend.position = "none")

plot_right <- ggplot(dataCollection, aes(aquVV)) + 
  geom_density(alpha=.5) + 
  coord_flip() + theme(axis.title.y = element_blank()) + #xlim(-30,0) +
  scale_y_continuous(breaks=c(0,0.2), labels=c("0","0.2"), limits=c(0,0.2)) +
  theme(legend.position = "none") 

grid.arrange(plot_top, empty, scatter, plot_right, ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4))
#ggsave(filename=path)
dev.off()

# world map accuracy plot

# calculate rmse by location
uniquePos <- unique(rslist[,6:7])
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

path <- paste(outdir, 'svr_rmse.png', sep="")

#dev.new(width=10, height=6, noRStudioGD=T)
png(filename = path, res=300, width = 2000, height = 1500)
gridPlt <- ggplot() + map + 
  geom_point(aes(x=lon, y=lat, colour=rmse), data=errors, size=0.8) + 
  scale_colour_gradientn(colours=r, name="RMSE", limit=c(0,0.3), na.value="red", breaks=c(0,0.1,0.2,0.3), labels=c("0","0.1","0.2","> 0.3"))
gridPlt + ggtitle("SVR prediction - RMSE\n")

dev.off()

path <- paste(outdir, 'svr_r.png', sep="")

#dev.new(width=10, height=6, noRStudioGD=T)
png(filename = path, res=300, width = 2000, height = 1500)
map <- borders("world", colour="gray50", fill="gray50")
#dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=r), data=errors, size=0.8) + 
  scale_colour_gradientn(colours=rev(r), name="R", labels=c("-1.0 ","-0.5","0.0","0.5","1.0"))
gridPlt + ggtitle("SVR prediction - R\n")

dev.off()
