#convert SMAP hdf files to GEOtiff
library(raster)
library(sp)
library(rhdf5)

SMAPtoGTiff <- function(inpath,outpath){
  
  lat <- as.vector(h5read(inpath, "/Sigma0_Data/cell_lat"))
  lon <- as.vector(h5read(inpath, "/Sigma0_Data/cell_lon"))
  vv_fore <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_vv_fore"))
  vv_aft <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_vv_aft"))
  hh_fore <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_hh_fore"))
  hh_aft <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_hh_aft"))
  vh_fore <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_xpol_fore"))
  vh_aft <- as.vector(h5read(inpath, "/Sigma0_Data/cell_sigma0_xpol_aft"))
  
  
  NAval <- -9999
  
  minlat <- 44
  minlon <- 4.5
  maxlat <- 50
  maxlon <- 17.5
  
  subsetInd <- which(lat >= minlat & lat <= maxlat & lon >= minlon & lon <= maxlon)
  lat <- lat[subsetInd]
  lon <- lon[subsetInd]
  vv_fore <- vv_fore[subsetInd]
  vv_aft <- vv_aft[subsetInd]
  hh_fore <- hh_fore[subsetInd]
  hh_aft <- hh_aft[subsetInd]
  vh_fore <- vh_fore[subsetInd]
  vh_aft <- vh_aft[subsetInd]
  
  vv_fore[vv_fore == -9999] <- NA
  vv_aft[vv_aft == -9999] <- NA
  hh_fore[hh_fore == -9999] <- NA
  hh_aft[hh_aft == -9999] <- NA
  vh_fore[vh_fore == -9999] <- NA
  vh_aft[vh_aft == -9999] <- NA
  
#   vv <- 10*log10((vv_fore+vv_aft)/2)
#   hh <- 10*log10((hh_fore+hh_aft)/2)
#   vh <- 10*log10((vh_fore+vh_aft)/2)
  vv <- 10*log10(vv_fore)
  hh <- 10*log10(hh_fore)
  vh <- 10*log10(vh_fore)
  
  rm(vv_fore,vv_aft,hh_fore,hh_aft,vh_fore,vh_aft)
  
  tmpl_raster <- raster(xmn=minlon, xmx=maxlon, ymn=minlat, ymx=maxlat, resolution=0.015)
  
  vv_raster <- rasterize(data.frame(x=lon,y=lat), tmpl_raster, field=vv, fun=mean, na.rm=TRUE)
  writeRaster(vv_raster, paste(outpath, "dB_vvfore.tif", sep=""), format="GTiff")
  hh_raster <- rasterize(data.frame(x=lon,y=lat), tmpl_raster, field=hh, fun=mean, na.rm=TRUE)
  writeRaster(hh_raster, paste(outpath, "dB_hhfore.tif", sep=""), format="GTiff")
  vh_raster <- rasterize(data.frame(x=lon,y=lat), tmpl_raster, field=vh, fun=mean, na.rm=TRUE)
  writeRaster(vh_raster, paste(outpath, "dB_vhfore.tif", sep=""), format="GTiff")
  
  
  H5close()
}