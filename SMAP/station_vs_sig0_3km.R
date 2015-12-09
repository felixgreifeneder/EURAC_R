#compare smap backscatter (3 km) with in-situ data. Based on data from Sab (only B3)

source("./SMAP/dailySMCavg.R")
library("RColorBrewer")
library("ggplot2")


#get SIG0
SIG0 <- read.table("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/B3SIG0_3km.txt",
                   header=T,
                   sep=",",
                   stringsAsFactors = F)

days <- as.Date(paste("2015",SIG0$mm,SIG0$dd,sep="-"))
SIG0 <- data.frame(date=days, SIG0)


#get in situ data
#b3 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB20.dat")
b3 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP20_scaled.dat")
b3$date <- as.Date(b3$date)

#create comparison array
VS <- data.frame(s0hh=SIG0$hh_linear, s0vv=SIG0$vv, ndist=SIG0$dist_nadir, smc=rep(-9999,length(days)))

for (i in 1:length(days)){
  
  if (any(b3$date == SIG0$date[i])){
    ind <- which(b3$date == SIG0$date[i])
    VS$s0hh[i] <- 10*log10(SIG0$hh_linear[i])
    VS$s0vv[i] <- 10*log10(SIG0$vv[i])
    VS$ndist[i] <- SIG0$dist_nadir[i]
    VS$smc[i] <- b3$mean_SM05[ind]
  }
  
}

plotVS <- data.frame(s0=c(VS$s0hh,VS$s0vv),
                     ndist=c(VS$ndist,VS$ndist),
                     smc=c(VS$smc,VS$smc),
                     pol=c(rep("HH",nrow(VS)),rep("VV",nrow(VS))))


val <- which(plotVS$ndist > 250)

p1 <- ggplot(data=plotVS[val,], aes(x=smc,y=s0 , colour=ndist)) + 
  geom_point()  +
  theme(aspect.ratio=1) +
  geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
  facet_grid(.~ pol) + ggtitle("in-situ SM vs SIG0 (3km), ndist>250\n")

dev.new()
p1

print("VV")
cor(VS$s0vv[VS$ndist>250], VS$smc[VS$ndist>250])^2
print("HH")
cor(VS$s0hh[VS$ndist>250], VS$smc[VS$ndist>250])^2