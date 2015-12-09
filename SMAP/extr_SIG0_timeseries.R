#extract SMAP SIG0 time-series for B1,B2, and B3

library(rhdf5)

extr_SMAP_SIG0_TS <- function(lat,lon){

  h5list <- list.files("D:/SMAP/L1C_SIG0/", pattern=".h5$", full.names=TRUE)
  SIG0 <- data.frame(vv_fore=as.numeric(), 
                     vv_aft=as.numeric(),
                    hh_fore=as.numeric(),
                    hh_aft=as.numeric(),
                    xpol_fore=as.numeric(),
                    xpol_aft=as.numeric(),
                    rflag=as.numeric(),
                    nadir_dist=as.numeric(),
                    elev_std=as.numeric(),
                    lia_fore=as.numeric(),
                    lia_aft=as.numeric(),
                    time=as.character(),
                    dist=as.numeric(),
                    lat=as.numeric(),
                    lon=as.numeric(),
                    stringsAsFactors = F)
  
  total <- length(h5list)
  pb <- txtProgressBar(min=0, max=total, style=3)

  for (fI in 1:length(h5list)){
  
      
    latlist <- h5read(h5list[fI], "/Sigma0_Data/cell_lat")
    lonlist <- h5read(h5list[fI], "/Sigma0_Data/cell_lon")
  
  
    arrayPos <- arrayInd(which.min(sqrt((latlist - lat)^2 + (lonlist - lon)^2)), .dim=c(nrow(latlist), ncol(latlist)))
  
    SIG0[fI,1] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_vv_fore", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,2] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_vv_aft", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,3] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_hh_fore", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,4] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_hh_aft", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,5] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_xpol_fore", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,6] <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_xpol_aft", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,7] <- h5read(h5list[fI], "/Sigma0_Data/cell_radar_mode_flag", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,8] <- h5read(h5list[fI], "/Crosstrack_Data/distance_from_nadir", index=list(x=arrayPos[1,1]))
    SIG0[fI,9] <- h5read(h5list[fI], "/Sigma0_Data/cell_altitude_std_dev", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,10] <- h5read(h5list[fI], "/Sigma0_Data/earth_incidence_mean_fore", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,11] <- h5read(h5list[fI], "/Sigma0_Data/earth_incidence_mean_aft", index=list(y=arrayPos[1,1], x=arrayPos[1,2]))
    SIG0[fI,12] <- h5read(h5list[fI], "/Spacecraft_Data/along_track_time_utc", index=list(x=arrayPos[1,2]))
    SIG0[fI,13] <- sqrt((latlist[arrayPos[1,1],arrayPos[1,2]] - lat)^2 + (lonlist[arrayPos[1,1],arrayPos[1,2]] - lon)^2)
    SIG0[fI,14] <- latlist[arrayPos[1,1],arrayPos[1,2]]
    SIG0[fI,15] <- lonlist[arrayPos[1,1],arrayPos[1,2]]
    
    if (sqrt((latlist[arrayPos[1,1],arrayPos[1,2]] - lat)^2 + (lonlist[arrayPos[1,1],arrayPos[1,2]] - lon)^2) < 1) {print(h5list[fI])}
    
    H5close()
    setTxtProgressBar(pb, fI)

  }
  return(SIG0)
  
}