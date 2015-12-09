#comparison between modeled soil moisture from ERA land and in situ smc (B1,B2,B3)
source("./SMAP/dailySMCavg.R")
source("./ERA/compute_daily_means.R")

#extracting the daily mean SMC from ERA, using the centre of B1, B2, B3 as a point of interest
eraDaily <- meanERA_SMC("D:/ERA/netcdf-atls03-a562cefde8a29a7288fa0b8b7f9413f7-oAQYlc.nc", c(46.6795983,10.58745841))

#extracting station data
b1 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB10.dat")
b2 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB15.dat")
b3 <- dailySMCavg("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB20.dat")
# b1 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP10_scaled.dat")
# b2 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP15_scaled.dat")
# b3 <- dailySMCavg_scaled("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSP20_scaled.dat")
b1$date <- as.Date(b1$date)
b2$date <- as.Date(b2$date)
b3$date <- as.Date(b3$date)

#format and collect data
dat <- data.frame(date=eraDaily$date, 
                  era=eraDaily$smc, 
                  b1=rep(-1, nrow(eraDaily)), 
                  b2=rep(-1, nrow(eraDaily)),
                  b3=rep(-1, nrow(eraDaily)))

for (i in 1:length(eraDaily$date)){
  
  if (any(b1$date == dat$date[i])) dat$b1[i] <- b1$mean_SM05[b1$date == dat$date[i]]
  if (any(b2$date == dat$date[i])) dat$b2[i] <- b2$mean_SM05[b2$date == dat$date[i]]
  if (any(b3$date == dat$date[i])) dat$b3[i] <- b3$mean_SM05[b3$date == dat$date[i]]
  
}

plot(dat$era, (dat$b1+dat$b2+dat$b3)/3)
cor(dat$era, (dat$b1+dat$b2+dat$b3)/3)
