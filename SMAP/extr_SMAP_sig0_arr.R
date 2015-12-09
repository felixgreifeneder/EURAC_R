#this function extracts SMAP sig0 data for a given geographical subset

library(rhdf5)

extr_SMAP_SIG0_arr <- function(path, minlat, minlon, maxlat, maxlon){
  
  latlist <- h5read(path, "/Sigma0_Data/cell_lat")
  lonlist <- h5read(path, "/Sigma0_Data/cell_lon")
  
  arrayPos <- which(latlist >= minlat & latlist <= maxlat &
                    lonlist >= minlon & lonlist <= maxlon, arr.ind=T)
  arrayIdx <- which(latlist >= minlat & latlist <= maxlat &
                      lonlist >= minlon & lonlist <= maxlon)
  
  SIG0 <- data.frame(x=numeric(length=nrow(arrayPos)),
                     y=numeric(length=nrow(arrayPos)),
                     s0vv_f=numeric(length=nrow(arrayPos)),
                     s0vv_a=numeric(length=nrow(arrayPos)),
                     s0hh_f=numeric(length=nrow(arrayPos)),
                     s0hh_a=numeric(length=nrow(arrayPos)),
                     s0xpol_f=numeric(length=nrow(arrayPos)),
                     s0xpol_a=numeric(length=nrow(arrayPos)),
                     date=character(length=nrow(arrayPos)),
                     stringsAsFactors = F)
  SIG0$x <- (h5read(path, "/Sigma0_Data/cell_lon"))[arrayIdx]
  SIG0$y <- (h5read(path, "/Sigma0_Data/cell_lat"))[arrayIdx]
  SIG0$s0vv_f <- (h5read(path, "/Sigma0_Data/cell_sigma0_vv_fore"))[arrayIdx]
  SIG0$s0vv_a <- (h5read(path, "/Sigma0_Data/cell_sigma0_vv_aft"))[arrayIdx]
  SIG0$s0hh_f <- (h5read(path, "/Sigma0_Data/cell_sigma0_hh_fore"))[arrayIdx]
  SIG0$s0hh_a <- (h5read(path, "/Sigma0_Data/cell_sigma0_hh_aft"))[arrayIdx]
  SIG0$s0xpol_f <- (h5read(path, "/Sigma0_Data/cell_sigma0_xpol_fore"))[arrayIdx]
  SIG0$s0xpol_a <- (h5read(path, "/Sigma0_Data/cell_sigma0_xpol_aft"))[arrayIdx]
  SIG0$date <- (h5read(path, "/Spacecraft_Data/along_track_time_utc"))[arrayPos[,2]]
  
  H5close()
  
  return(SIG0)
  
  
}