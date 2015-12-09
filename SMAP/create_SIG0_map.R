library("rhdf5")
library("RColorBrewer")
library("ggplot2")
library("ggmap")
library("scales")

smapL1filelist <- list.files(path="X:/Workspaces/GrF/01_Data/SMAP/L1C_SIG0/", pattern="*.h5$", full.names=T)

for (fInd in 1:length(smapL1filelist)){

  lon <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_lon")
  lat <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_lat")
  sig0hh_aft <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_sigma0_hh_aft")
  sig0hh_fore <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_sigma0_hh_fore")
  sig0vv_aft <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_sigma0_vv_aft")
  sig0vv_fore <- h5read(smapL1filelist[fInd], "/Sigma0_Data/cell_sigma0_vv_fore")
  
  doa <- substr(smapL1filelist[fInd],59,66)
  fname <- strsplit(basename(smapL1filelist[fInd]), ".", fixed=T)
  fname <- fname[[1]][1]

  NAval <- -9999

  minlat <- 44
  minlon <- 4.5
  maxlat <- 50
  maxlon <- 17.5

  val <- which(sig0vv_fore != NAval & (lat >= minlat & lat <= maxlat) & (lon >= minlon & lon <= maxlon))
  if (length(val) == 0) {next}
  tmpDf <- data.frame(lat=lat[val], lon=lon[val], sig0=sig0vv_fore[val], sig0stretch=sig0vv_fore[val])

  rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r <- rf(32)
  
  #create colour stretch
  val_str <- 1.5*sd(tmpDf$sig0)
  vInd <- which(tmpDf$sig0>=(mean(tmpDf$sig0)-val_str) & tmpDf$sig0<=(mean(tmpDf$sig0)+val_str))
  #vInd <- which(tmpDf$bt>=220)
  lInd <- which(tmpDf$sig0<(mean(tmpDf$sig0)-val_str))
  #lInd <- which(tmpDf$bt<220)
  hInd <- which(tmpDf$sig0>(mean(tmpDf$sig0)+val_str))
  #hInd <- which(tmpDf$bt>max(tmpDf$bt))
  tmpDf$sig0stretch[vInd] <- rescale(tmpDf$sig0[vInd])
  tmpDf$sig0stretch[lInd] <- 0
  tmpDf$sig0stretch[hInd] <- 1
  
  #sort
#   sInd <- sort(tmpDf$bt, index.return=T)
#   sInd <- sInd$ix
#   tmpDf <- tmpDf[sInd,]
  

  map <- get_map(location=c(minlon,minlat,maxlon,maxlat), maptype="terrain", source="google")
  dev.new(width=9, height=8, noRStudioGD=T)
  gridPlt <- ggmap(map) + 
    geom_point(aes(x=lon, y=lat, colour=sig0stretch), data=tmpDf, size=0.4) + 
    scale_colour_gradientn(colours=r, name="SIG0")
  gridPlt + ggtitle(paste("SMAP-",doa))
  ggsave(file = paste("X:/Workspaces/GrF/01_Data/SMAP/L1C_SIG0/qlooks/", fname, ".png", sep=""))
  dev.off()
  H5close()
}
