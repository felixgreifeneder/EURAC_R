bdown=function(url, file){
  library('RCurl')
  f = CFILE(file, mode="wb")
  a = curlPerform(url = url, writedata = f@ref, noprogress=FALSE)
  close(f)
  return(a)
}


#check files ... find files in extent
library("raster")
library("rgdal")
library("rhdf5")

filelist <- read.table("X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/filelist_D_052913.txt", as.is=T)
proj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ")
point_values <- data.frame()

for (name in filelist$V1){
  destF <- paste("X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/", basename(name), sep="")
  if (file.exists(destF) == F){
    print(paste("Downloading", name))
    #download.file(url=name, destfile=destF, mode="wb")
    ret = bdown(name, destF)
    print(paste("Errors ", ret))
  }
  #latmin, latmax, lonmin, lonmax
  print("Checking extent ...")
  bbox <- c(46.6455, 46.7795, 10.5601, 10.7371)
  
  h5lat <- h5read(destF, "Sigma0_Data/cell_lat")
  h5lon <- h5read(destF, "Sigma0_Data/cell_lon")
  val <- which(h5lat != -9999 & h5lat >= bbox[1] & h5lat <= bbox[2] & 
               h5lon != -9999 & h5lon >= bbox[3] & h5lon <= bbox[4])
  
  if (length(val) > 0){
    print("Inside!")
    
    #extract data and resample
    h5lat_subset <- h5lat[val]
    h5lon_subset <- h5lon[val]
    if (exists("pointValshh") == F){
      pointValshh <- data.frame(x=h5lon_subset, y=h5lat_subset)
      pointValsvv <- data.frame(x=h5lon_subset, y=h5lat_subset)
      pointValshv <- data.frame(x=h5lon_subset, y=h5lat_subset)
    }
    xdim <- ceiling((max(h5lon_subset)-min(h5lon_subset)) / 0.013)
    xmin <- min(h5lon_subset)
    xmax <- xmin + (xdim*0.013)
    ydim <- ceiling((max(h5lat_subset)-min(h5lat_subset)) / 0.013)
    ymin <- min(h5lat_subset)
    ymax <- ymin + (ydim*0.013)
    
    init_mat <- matrix(rep(-Inf, (xdim*ydim)), ydim, xdim)
    init_raster <- raster(init_mat, xmn=xmin, xmx=xmax, ymn=ymin, ymx=ymax, crs=proj)
    NAvalue(init_raster) <- -Inf
    hh_aft <- init_raster
    h5hh <- h5read(destF, "Sigma0_Data/cell_sigma0_hh_aft")
    h5hh_subset <- h5hh[val]
    pointValshh[[substr(basename(name), 27, 34)]] <- h5hh_subset
    val2 <- which(h5hh_subset != -9999)
    xy <- cbind(h5lon_subset[val2], h5lat_subset[val2])
    hh_aft <- rasterize(xy, hh_aft, h5hh_subset[val2], background=NA, fun=mean)
    rm(h5hh, h5hh_subset)
    
    vv_aft <- init_raster
    h5vv <- h5read(destF, "Sigma0_Data/cell_sigma0_vv_aft")
    h5vv_subset <- h5vv[val]
    pointValsvv[[substr(basename(name), 27, 34)]] <- h5vv_subset
    val2 <- which(h5vv_subset != -9999)
    xy <- cbind(h5lon_subset[val2], h5lat_subset[val2])
    vv_aft <- rasterize(xy, vv_aft, h5vv_subset[val2], background=NA, fun=mean)
    rm(h5vv, h5vv_subset)
    
    hv_aft <- init_raster
    h5hv <- h5read(destF, "Sigma0_Data/cell_sigma0_xpol_aft")
    h5hv_subset <- h5hv[val]
    pointValshv[[substr(basename(name), 27, 34)]] <- h5hv_subset
    val2 <- which(h5hv_subset != -9999)
    xy <- cbind(h5lon_subset[val2], h5lat_subset[val2])
    hv_aft <- rasterize(xy, hv_aft, h5hv_subset[val2], background=NA, fun=mean)
    rm(h5hv, h5hv_subset)
    
    writeRaster(hh_aft, paste("X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/", substr(basename(name), 1, 52), "_hh.tif", sep=""), "GTiff", overwrite=TRUE)
    writeRaster(vv_aft, paste("X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/", substr(basename(name), 1, 52), "_vv.tif", sep=""), "GTiff", overwrite=TRUE)
    writeRaster(hv_aft, paste("X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/", substr(basename(name), 1, 52), "_hv.tif", sep=""), "GTiff", overwrite=TRUE)
    H5close()
    file.remove(destF)
  } else {
    print("Outside")
    file.remove(destF)
  }
}

#convert to SpatialPointsDataFrame
sphh <- SpatialPointsDataFrame(pointValshh[c("x", "y")], proj4string=proj, data=pointValshh[,c(3:22)])
spvv <- SpatialPointsDataFrame(pointValsvv[c("x", "y")], proj4string=proj, data=pointValsvv[,c(3:22)])
sphv <- SpatialPointsDataFrame(pointValshh[c("x", "y")], proj4string=proj, data=pointValshv[,c(3:22)])

writeOGR(sphh, "X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/SMAP_grid_SIG0_HH.shp", layer="SIG0_hh", driver="ESRI Shapefile")
writeOGR(spvv, "X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/SMAP_grid_SIG0_VV.shp", layer="SIG0_vv", driver="ESRI Shapefile")
writeOGR(sphv, "X:/Workspaces/GrF/01_Data/SMAP/sig0_simulations/SMAP_grid_SIG0_HV.shp", layer="SIG0_hv", driver="ESRI Shapefile")

