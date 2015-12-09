#extract SMAP SIG0 time-series for B1,B2, and B3

library(rhdf5)
library(MASS)

extr_SMAP_SIG0_TS_bilin <- function(lat,lon){

  h5list <- list.files("D:/SMAP/L1C_SIG0/", pattern=".h5$", full.names=TRUE)
  nSamples <- length(h5list)
  SIG0 <- data.frame(vv_fore=rep(-9999, nSamples), 
                     vv_aft=rep(-9999, nSamples),
                     hh_fore=rep(-9999, nSamples),
                     hh_aft=rep(-9999, nSamples),
                     xpol_fore=rep(-9999, nSamples),
                     xpol_aft=rep(-9999, nSamples),
                     rflag=rep(-9999, nSamples),
                     nadir_dist=rep(-9999, nSamples),
                     elev_std=rep(-9999, nSamples),
                     lia_fore=rep(-9999, nSamples),
                     lia_aft=rep(-9999, nSamples),
                     time=rep("2015-01-01", nSamples),
                     ease_x=rep(-9999, nSamples),
                     ease_y=rep(-9999, nSamples),
                     stringsAsFactors = F)
  
  pb <- txtProgressBar(min=0, max=nSamples, style=3)

  for (fI in 1:nSamples){
  
    latlist <- h5read(h5list[fI], "/Sigma0_Data/cell_lat")
    lonlist <- h5read(h5list[fI], "/Sigma0_Data/cell_lon")
    rflaglist <- h5read(h5list[fI], "/Sigma0_Data/cell_radar_mode_flag")
    
    #calculate distances to all pixels within the grid and select the closest 4 
    #for bilinear interpolation
    distances <- sqrt((latlist - lat)^2 + (lonlist - lon)^2)
    sort_dist <- sort(distances, index.return=T)
    sort_ind <- sort_dist[["ix"]]
    #print(distances[sort_ind[1]])
    # if distance > 1km, skip to the next file
    if (distances[sort_ind[1]] >= 0.01) {
      H5close()
      setTxtProgressBar(pb, fI)
      next
    }
    
    best4ind <- sort_ind[1:4]
    best4XYind <- arrayInd(sort_ind[1:4], .dim = c(nrow(latlist), ncol(latlist)))
    best4lat <- latlist[best4ind]
    best4lon <- lonlist[best4ind]
    rflag <- rflaglist[best4ind[1]]
    
    #determin the weight vector
    w <- matrix(distances[best4ind]/sum(distances[best4ind]), ncol=4, nrow=1)
    w <- matrix(w[1,4:1], ncol=4, nrow = 1)
    
    #read each variable and apply bilinear interpolation
    #vv_fore
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_vv_fore")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$vv_fore[fI] <- w %*% b
    } else {
      SIG0$vv_fore[fI] <- -9999
    }
    #vv_aft
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_vv_aft")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$vv_aft[fI] <- w %*% b
    } else {
      SIG0$vv_aft[fI] <- -9999
    }
    #hh_fore
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_hh_fore")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$hh_fore[fI] <- w %*% b
    } else {
      SIG0$hh_fore[fI] <- -9999
    }
    #hh_aft
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_hh_aft")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$hh_aft[fI] <- w %*% b
    } else {
      SIG0$hh_aft[fI] <- -9999
    }
    #xpol_fore
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_xpol_fore")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$xpol_fore[fI] <- w %*% b
    } else {
      SIG0$xpol_fore[fI] <- -9999
    }
    #xpol_aft
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_sigma0_xpol_aft")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$xpol_aft[fI] <- w %*% b
    } else {
      SIG0$xpol_aft[fI] <- -9999
    }
    #nadir_dist
    tmp <- h5read(h5list[fI], "/Crosstrack_Data/distance_from_nadir")
    best4val <- tmp[best4XYind[,1]]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$nadir_dist[fI] <- w %*% b
    } else {
      SIG0$nadir_dist[fI] <- -9999
    }
    #elev_std
    tmp <- h5read(h5list[fI], "/Sigma0_Data/cell_altitude_std_dev")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$elev_std[fI] <- w %*% b
    } else {
      SIG0$elev_std[fI] <- -9999
    }
    #lia_fore
    tmp <- h5read(h5list[fI], "/Sigma0_Data/earth_incidence_mean_fore")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$lia_fore[fI] <- w %*% b
    } else {
      SIG0$lia_fore[fI] <- -9999
    }
    #lia_aft
    tmp <- h5read(h5list[fI], "/Sigma0_Data/earth_incidence_mean_aft")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >= 16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$lia_aft[fI] <- w %*% b
    } else {
      SIG0$lia_aft[fI] <- -9999
    }
    #time
    tmp <- h5read(h5list[fI], "Spacecraft_Data/along_track_time_utc")
    best4val <- tmp[best4XYind[,2]]
    SIG0$time[fI] <- best4val[1]
    #ease grid row index
    tmp <- h5read(h5list[fI],"/Sigma0_Data/cylindrical_grid_row_index")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >=16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$ease_y[fI] <- w %*% b
    } else {
      SIG0$ease_y[fI] <- -9999    
    }
    #ease grid col index
    tmp <- h5read(h5list[fI],"/Sigma0_Data/cylindrical_grid_column_index")
    best4val <- tmp[best4ind]
    if (all(best4val != -9999) & (rflag < 8 | rflag >=16)){
      b <- matrix(data=best4val, ncol=1, nrow=4)
      SIG0$ease_x[fI] <- w %*% b
    } else {
      SIG0$ease_x[fI] <- -9999    
    }
  
    rm(tmp)
 
    H5close()
    setTxtProgressBar(pb, fI)

  }
  
  SIG0 <- SIG0[SIG0$vv_fore != -9999,]
  return(SIG0)
  
}