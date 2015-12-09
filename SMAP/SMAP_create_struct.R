#resample and collect all information in one data frame

library(rhdf5)

#create lat lon grid

minlat <- 44
minlon <- 4.5
maxlat <- 50 
maxlon <- 17.5

latlist <- seq(minlat, maxlat, 0.015)
lonlist <- seq(minlon, maxlon, 0.015)
nelements <- length(latlist)*length(lonlist)

grid <- data.frame(locID=numeric(length=nelements), 
                   lat=numeric(length=nelements), 
                   lon=numeric(length=nelements))

pb <- txtProgressBar(min=0, max=nelements, style=3)
cntr <- 1
for (i in 1:length(latlist)){
  for (j in 1:length(lonlist)){
    grid$locID[cntr] <- cntr
    grid$lat[cntr] <- latlist[i]
    grid$lon[cntr] <- lonlist[j]
    cntr <- cntr + 1
  }
  setTxtProgressBar(pb, cntr)
}

save(grid, file = "./grid.dat")

flist <- list.files("D:/SMAP/L1C_SIG0", pattern=".h5$", full.names=TRUE)
total <- length(flist)
SMAPdata <- data.frame(day=character(length=total*nrow(grid)),
                       locID=numeric(length=total*nrow(grid)),
                       sig0vv_f=numeric(length=total*nrow(grid)),
                       sig0vv_a=numeric(length=total*nrow(grid)),
                       sig0hh_f=numeric(length=total*nrow(grid)),
                       sig0hh_a=numeric(length=total*nrow(grid)),
                       sig0xpol_f=numeric(length=total*nrow(grid)),
                       sig0xpol_a=numeric(length=total*nrow(grid)), stringsAsFactors = F)


pb <- txtProgressBar(min=0, max=total, style=3)
cntr <- 1
for (fID in 1:length(flist)){
  
  SMAPfile <- extr_SMAP_SIG0_arr(flist[fID], minlat=44,minlon=4.5,maxlat=50,maxlon=17.5)
  for (i in 1:nrow(SMAPfile)){
    dist <- sqrt((SMAPfile$x[i]-grid$lon)^2 + (SMAPfile$y[i]-grid$lat)^2)
    if (min(dist) > 0.015 | any(is.finite(dist)) == F) next
    if (SMAPfile$s0vv_f[i] == -9999) next
    minInd <- which.min(dist)
    SMAPdata$day[cntr] <- as.Date(SMAPfile$date[i])
    SMAPdata$locID[cntr] <- minInd
    SMAPdata$sig0vv_f[cntr] <- SMAPfile$s0vv_f[i]
    SMAPdata$sig0vv_a[cntr] <- SMAPfile$s0vv_a[i]
    SMAPdata$sig0hh_f[cntr] <- SMAPfile$s0hh_f[i]
    SMAPdata$sig0hh_a[cntr] <- SMAPfile$s0hh_a[i]
    SMAPdata$sig0xpol_f[cntr] <- SMAPfile$s0xpol_f[i]
    SMAPdata$sig0xpol_a[cntr] <- SMAPfile$s0xpol_a[i]
    cntr <- cntr+1
    
  }
  
  setTxtProgressBar(pb, fID)
  H5close()
  
}