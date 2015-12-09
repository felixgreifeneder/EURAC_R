#derive SMC map

library(raster)
source("./Aquarius/parallel_predictions.R")

workdir = "C:/Users/FGreifeneder/Documents/work/"

#load and prepared datasets
rs2path <- "X:/Workspaces/GrF/01_Data/RADARSAT2/2015/RS2_OK64176_PK603574_DK535224_S3_20150709_170657_HH_HV_SLC.tif"
file.copy(rs2path, paste(workdir, basename(rs2path), sep=""))
rs2path <- paste(workdir, basename(rs2path), sep="")
rs2hh <- raster(rs2path, band=1)
rs2hv <- raster(rs2path, band=2)
rs2lia <- raster(rs2path, band=3)

#load dem
dempath <- "X:/Workspaces/GrF/01_Data/ANCILLARY/DEM/dtm_10_st_sub_west_latlon.tif"
file.copy(dempath, paste(workdir, basename(dempath), sep=""))
dempath <- paste(workdir, basename(dempath), sep="")
dem <- raster(dempath)
slopepath <- "X:/Workspaces/GrF/01_Data/ANCILLARY/DEM/dtm_10_st_sub_west_latlon_slope.tif"
file.copy(slopepath, paste(workdir, basename(slopepath), sep=""))
slopepath <- paste(workdir, basename(slopepath), sep="") 
slope <- raster(slopepath)
aspectpath <- "X:/Workspaces/GrF/01_Data/ANCILLARY/DEM/dtm_10_st_sub_west_latlon_aspect.tif"
file.copy(aspectpath, paste(workdir, basename(aspectpath), sep=""))
aspectpath <- paste(workdir, basename(aspectpath), sep="")
aspect <- raster(aspectpath)

#crop datasets
ROIext <- extent(10.51, 10.70, 46.634, 46.808)
rs2hh <- crop(rs2hh, ROIext)
rs2hv <- crop(rs2hv, ROIext)
rs2lia <- crop(rs2lia, ROIext)
dem <- crop(dem, ROIext)
slope <- crop(slope, ROIext)
aspect <- crop(aspect, ROIext)


#--------------------------------------------------------------------------------------------------
#extract data

df <- data.frame(hh=numeric(),
                 hv=numeric(),
                 lia=numeric(),
                 height=numeric(),
                 slope=numeric(),
                 aspect=numeric(),
                 row=numeric(),
                 col=numeric())

cntr <- 1

pb <- txtProgressBar(min = 0, max = nrow(rs2hh)*ncol(rs2hh))

for (rind in 1:nrow(rs2hh)){
  for (cind in 1:ncol(rs2hh)){
    
    setTxtProgressBar(pb, value=cntr)
    lon <- xFromCol(rs2hh, col=cind)
    lat <- yFromRow(rs2hh, row=rind)
    lonlat <- matrix(c(lon,lat),1,2)
    
    hh <- rs2hh[rind,cind]
    hv <- rs2hv[rind,cind]
    lia <- rs2lia[rind,cind]
    h <- extract(dem, lonlat)
    s <- extract(slope, lonlat)
    a <- extract(aspect, lonlat)
    r <- rind
    c <- cind
    
    tmp <- data.frame(hh=hh,
                      hv=hv,
                      lia=lia,
                      height=h,
                      slope=s,
                      aspect=a,
                      row=r,
                      col=c)
    df[cntr,] <- tmp
    cntr <- cntr + 1
    
    
  }
}
close(pb)

#--------------------------------------------------------------------------------------------------
#SVR
#--------------------------------------------------------------------------------------------------

load("./Radarsat2Mazia/svrModel.dat")
tunedModel <- SVRtuning$best.model

#initiate output map
smcRaster <- raster(rs2hh)

#preditc smc
SMCpredicted <- parallel_predictions(tunedModel, df)

#fill values into new map
for (i in 1:nrow(df)){
  valid <- df$lia[i] > 10 & df$lia[i] < 50
  
  if (valid == T){
    SMCRaster[df$row[i], df$col[i]] <- SMCpredicted[i]
  } else {
    SMCRaster[df$row[i], df$col[i]] <- -1
  }
}

writeRaster(SMCRaster, paste(workdir, "SMCmap.tif", seP=""))
