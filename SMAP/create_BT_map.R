library("rhdf5")
library("RColorBrewer")
library("ggplot2")
library("ggmap")
library("scales")

smapL1filelist <- list.files(path="X:/Workspaces/GrF/01_Data/SMAP/L1C_TB/", pattern="*.h5", full.names=T)

for (fInd in 1:length(smapL1filelist)){

  lon <- h5read(smapL1filelist[fInd], "/Global_Projection/cell_lon")
  lat <- h5read(smapL1filelist[fInd], "/Global_Projection/cell_lat")
  bt <- h5read(smapL1filelist[fInd], "/Global_Projection/cell_tb_v_fore")
  
  doa <- substr(smapL1filelist[fInd],59,66)
  fname <- strsplit(basename(smapL1filelist[fInd]), ".", fixed=T)
  fname <- fname[[1]][1]

  NAval <- -9999

  minlat <- 44
  minlon <- 4.5
  maxlat <- 50
  maxlon <- 17.5

  val <- which(bt != NAval & (lat >= minlat & lat <= maxlat) & (lon >= minlon & lon <= maxlon))
  if (length(val) == 0) {next}
  tmpDf <- data.frame(lat=lat[val], lon=lon[val], bt=bt[val], btstretch=rep(0, length(bt[val])))

  rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
  r <- rf(32)
  
  #create colour stretch
  val_str <- sd(tmpDf$bt)
  #vInd <- which(tmpDf$bt>=(mean(tmpDf$bt)-val_str) & tmpDf$bt<=(mean(tmpDf$bt)+val_str))
  vInd <- which(tmpDf$bt>=220)
  #lInd <- which(tmpDf$bt<(mean(tmpDf$bt)-val_str))
  lInd <- which(tmpDf$bt<220)
  #hInd <- which(tmpDf$bt>(mean(tmpDf$bt)+val_str))
  #hInd <- which(tmpDf$bt>max(tmpDf$bt))
  tmpDf$btstretch[vInd] <- rescale(tmpDf$bt[vInd])
  tmpDf$btstretch[lInd] <- 0
  #tmpDf$btstretch[hInd] <- 1
  
  #sort
  sInd <- sort(tmpDf$bt, index.return=T)
  sInd <- sInd$ix
  tmpDf <- tmpDf[sInd,]
  

  map <- get_map(location=c(minlon,minlat,maxlon,maxlat), maptype="terrain", source="google")
  dev.new(width=9, height=8, noRStudioGD=T)
  gridPlt <- ggmap(map) + 
    geom_point(aes(x=lon, y=lat, colour=btstretch), data=tmpDf, size=8) + 
    scale_colour_gradientn(colours=r, name="K")
  gridPlt + ggtitle(paste("SMAP-",doa))
  ggsave(file = paste("X:/Workspaces/GrF/01_Data/SMAP/L1C_TB/qlooks/", fname, ".png", sep=""))
  dev.off()
}
