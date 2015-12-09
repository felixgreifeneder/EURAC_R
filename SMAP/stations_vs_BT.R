#compare stations and radiometer BT

library("rhdf5")
library("ggplot2")
library("gtable")
library("grid")

#create the BT time-series

radio_files <- list.files(path="X:/Workspaces/GrF/01_Data/SMAP/L1C_TB/", pattern="*.h5", full.names=T)

#center between B1, B2, B3 as a target location
cLon <- 10.587
cLat <- 46.68

btTS <- data.frame()

#find the closest SMAP grid point in relation to the staion locations
#and store brightness temperature time-series
for (radio_path in radio_files){
  
  lon <- h5read(radio_path, "/Global_Projection/cell_lon")
  lat <- h5read(radio_path, "/Global_Projection/cell_lat")
  btV_f <- h5read(radio_path, "/Global_Projection/cell_tb_v_fore")
  btV_a <- h5read(radio_path, "/Global_Projection/cell_tb_v_aft")
  btH_f <- h5read(radio_path, "/Global_Projection/cell_tb_h_fore")
  btH_a <- h5read(radio_path, "/Global_Projection/cell_tb_h_aft")
  time <- (h5read(radio_path, "/Global_Projection/cell_tb_time_seconds_aft") +
          h5read(radio_path, "/Global_Projection/cell_tb_time_seconds_fore")) / 2
  
  dist <- sqrt((lon-cLon)^2+(lat-cLat^2))
  minDist <- which.min(dist)
  
  btTS <- rbind(btTS, data.frame(time=time[minDist],
                                 lon=lon[minDist], 
                                 lat=lat[minDist], 
                                 btV_a=btV_a[minDist],
                                 btV_f=btV_f[minDist],
                                 btH_a=btH_a[minDist],
                                 btH_f=btH_f[minDist]))
  
}

#convert J2000 seconds to posix object
posixOrigin <- as.POSIXct("2000-01-01-11:58:56", format="%F-%H:%M:%S", tz="UTC")
btTS$time <- as.POSIXct(btTS$time, origin=posixOrigin, tz="UTC")

stTS <- data.frame()


load("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB10.dat")
b1 <- combinedTable
load("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB15.dat")
b2 <- combinedTable
load("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/FullTSB20.dat")
b3 <- combinedTable

stTS <- data.frame(station=rep("b1", nrow(b1)), 
                     time=as.POSIXct(paste(b1$Yr,"-",b1$Mo,"-",b1$Dy,"-",b1$Hr,":",b1$Min,sep=""), format="%F-%H:%M", tz="CET"), b1)
stTS <- rbind(b1b2b3, data.frame(station=rep("b2", nrow(b2)),
                                   time=as.POSIXct(paste(b2$Yr,"-",b2$Mo,"-",b2$Dy,"-",b2$Hr,":",b2$Min,sep=""), format="%F-%H:%M", tz="CET"), b2))
stTS <- rbind(b1b2b3, data.frame(station=rep("b3", nrow(b3)),
                                   time=as.POSIXct(paste(b3$Yr,"-",b3$Mo,"-",b3$Dy,"-",b3$Hr,":",b3$Min,sep=""), format="%F-%H:%M", tz="CET"), b3))
rm(b1,b2,b3)

valid <- which(stTS$SM05_1 != -7777 & stTS$SM05_2 != -7777)# & stTS$Yr==2015 & stTS$Mo>=7)
stTS <- stTS[valid,]
posixFirstday <- as.POSIXct("2015-07-01-00:00:00", format="%F-%H:%M:%S", tz="UTC")
validBT <- which(btTS$time > posixFirstday)
btTS <- btTS[validBT,]

p1 <- ggplot(data=stTS, aes(x=time,y=((SM05_1+SM05_2)/2),group=station,colour=station)) + geom_line() 
p2 <- ggplot(data=btTS, aes(x=time,y=btV_a)) + geom_point()

dev.new()
p1
dev.new()
p2




