# this script reads 

library(raster)
library(ncdf4)

createRasterforExtent <- function(extent, res){
  
  minx <- extent[1]
  miny <- extent[2]
  maxx <- extent[3]
  maxy <- extent[4]
  
  outraster <- raster(xmn=minx, xmx=maxx, ymn=miny, ymx=maxy, resolution=res)
  
  return(outraster)
  
}


findNN <- function(latlist, lonlist, gpilist, px){
  
  min <- 9999
  gpiid <- NA
  
  for (i in 1:length(gpilist)){
    
    dx <- px[1] - lonlist[i]
    dy <- px[2] - latlist[i]
    
    dist <- sqrt(dx^2+dy^2)
    
    if (dist < min){ 
      min <- dist 
      gpiid <- gpilist[i]
    }
    
  }
  
  return(gpiid)
  
}


gpiPList <- function(inraster, warp5grid_path, extent){
  
  # load ncdef
  warp5grid_def <- nc_open(warp5grid_path)
  
  gpilist <- ncvar_get(warp5grid_def, varid = 'gpi')
  latlist <- ncvar_get(warp5grid_def, varid = 'lat')
  lonlist <- ncvar_get(warp5grid_def, varid = 'lon')
  
  # crop grid to extent
  validgpis <- which(lonlist >= extent[1]-1 & lonlist <= extent[3]+1 & latlist >= extent[2]-1 & latlist <= extent[4]+1)
  gpilist <- gpilist[validgpis]
  latlist <- latlist[validgpis]
  lonlist <- lonlist[validgpis]
  
  # initialising list relating each pixel of the raster
  # with a gpi
  gridPlist <- c()
  
  # iterating through each pixel of the raster
  for (ix in c(1:ncol(inraster))){
    for (iy in c(1:nrow(inraster))){
      
      ilon <- xFromCol(inraster, col=ix)
      ilat <- yFromRow(inraster, row=iy)
      
      gpiid <- findNN(latlist, lonlist, gpilist, c(ilon,ilat))
      gridPlist <- c(gridPlist, gpiid)
      
    }
  }
  
  return(gridPlist)
  
}

createImage <- function(rastertempl, tsfile, date_oi, gpilist){
  
  ts_nc <- nc_open(tsfile)
  # load location ids
  loc_id <- ncvar_get(ts_nc, varid = 'location_id')
  # load row size
  row_size_all <- ncvar_get(ts_nc, varid = 'row_size')
  # load sm
  sm_all <- ncvar_get(ts_nc, varid = 'sm')
  # load dates
  dates_all <- ncvar_get(ts_nc, varid = 'time')
  # load surface state
  frozen_all <- ncvar_get(ts_nc, varid = 'ssf')
  # wetland/water
  #wet_all <- ncvar_get(ts_nc, varid = 'advf_wetland')
  
  outraster <- rastertempl
  NAvalue(outraster) <- -1
  
  #iterate through each pixel
  pxcounter <- 1
  
  for (ix in 1:ncol(outraster)){
    for (iy in 1:nrow(outraster)){
      
      tsid <- which(loc_id == gpilist[pxcounter])
      if (length(tsid) == 0){
        sm_out <- -1
      } else {
        row_size <- row_size_all[tsid]
        row_start <- sum(row_size_all[1:tsid-1])+1
        row_end <- row_start + row_size-1
      
        # extract current time-series
        dates_ts <- trunc(dates_all[row_start:row_end])
        sm_ts <- sm_all[row_start:row_end]
        frozen_ts <- frozen_all[row_start:row_end]
        #wet_ts <- dates_all[row_start:row_end]
      
        # extract measurement at day of interest
        did <- which(dates_ts == date_oi)
        if (length(did) == 0){
          sm_out <- -1
        } else {
          sm_single <- sm_ts[did]
          if (length(sm_single) > 1) { sm_single <- sm_single[1] }
          if (is.na(sm_single) == TRUE) { sm_single <- 127 }
          frozen_single <- frozen_ts[did]
          if (length(frozen_single) > 1) {frozen_single <- frozen_single[1]}
          if (is.na(frozen_single == TRUE)) {frozen_single <- 0 }
          #wet_single <- wet_ts[did]
      
          # check if sm is in valid range or any of the flags are set
          if (sm_single == 127 | frozen_single != 1){
            sm_out <- -1
          } else {
            sm_out <- sm_single
          }
        }
      }
      
      outraster[iy,ix] <- sm_out
      
      pxcounter <- pxcounter + 1
      
    }
  }
  
  return(outraster)
  
}


extent <- c(20.4, 40.5, 21.5, 41.5)
resolution <- 0.01
outpath <- 'T:\\ECOPOTENTIAL\\MAC-AL - LAKE OHRID\\SM_ASCAT_GTIFS\\'

# initialise raster
rastertempl <- createRasterforExtent(extent, resolution)
grd_pxl_list <- gpiPList(rastertempl, "T:\\ECOPOTENTIAL\\MAC-AL - LAKE OHRID\\TUW_WARP5_grid_info_2_1.nc", extent)

#cycle throug dates
startdate <- as.Date("2013-06-25")
enddate <- as.Date("2015-12-31")

for (idate in startdate:enddate){
  print(idate+25567)
  print(as.character(as.Date(idate, origin='1970-01-01'), format="%Y%m%d"))
  traster <- createImage(rastertempl, "T:\\ECOPOTENTIAL\\MAC-AL - LAKE OHRID\\SM_ASCAT_TS12.5_DR2016\\1466.nc", idate+25567, grd_pxl_list)
  if (length(unique(traster)) > 4){
    filepath <- paste(outpath, 'SM_ASCAT_H109_', as.character(as.Date(idate, origin='1970-01-01'), format="%Y%m%d"), '.tif', sep="")
    writeRaster(traster, filename = filepath, format = "GTiff")
  }
}
