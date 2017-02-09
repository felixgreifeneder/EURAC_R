#script for the production of smc maps

library(raster)
library(gbm)
library(caret)
library(e1071)
library(ggplot2)
library(foreach)
library(doSNOW)

#working dir
work <- 'C:/Users/FGreifeneder/Documents/work'

#-------------------------------------------------------------
#Initiate required datasets
#-------------------------------------------------------------

#SAR dataset
rs2hhpath <- "X:/ProjectData/HiResAlp/04_Data/SM SVR Argentina/Images/RS2-SLC-S4-DES-09-Aug-2015_HH_ML_Spk_TC_band1.tif"
rs2hvpath <- "X:/ProjectData/HiResAlp/04_Data/SM SVR Argentina/Images/RS2-SLC-S4-DES-09-Aug-2015_HV_ML_Spk_TC_band1.tif"
rs2liapath <- "X:/ProjectData/HiResAlp/04_Data/SM SVR Argentina/LIA/09-Ago-2015/incidenceAngle_HH.img"
rs2liahdr <- paste(dirname(rs2liapath), "/", substr(basename(rs2liapath), 1, nchar(basename(rs2liapath))-4), ".hdr", sep="")

file.copy(rs2hhpath, paste(work, basename(rs2hhpath), sep="/"))
file.copy(rs2hvpath, paste(work, basename(rs2hvpath), sep="/"))
file.copy(rs2liapath, paste(work, basename(rs2liapath), sep="/"))
file.copy(rs2liahdr, paste(work, basename(rs2liahdr), sep="/"))

rs2hhpath <- paste(work, basename(rs2hhpath), sep="/")
rs2hvpath <- paste(work, basename(rs2hvpath), sep="/")
rs2liapath <- paste(work, basename(rs2liapath), sep="/")

rs2hh <- raster(rs2hhpath)
rs2hv <- raster(rs2hvpath)
rs2lia <- raster(rs2liapath)
##crop
ext <- extent(matrix(c(-62.52,-33.05,-62.47,-33.0), ncol = 2, nrow=2))
rs2hh <- crop(rs2hh, ext)
rs2hv <- resample(rs2hv, rs2hh)
rs2lia <- resample(rs2lia, rs2hh)

#NDVI
NDVIpath <- "C:/Users/FGreifeneder/Documents/tmp_proc/RS2_argentina/NDVI/MOD13Q1_A2015209_h12v12_005_2015226081623_250m_16_days_NDVI_3154f5a0.tif"
NDVI <- raster(NDVIpath)
NDVIres <- resample(NDVI, rs2hh, method="bilinear")
remove(NDVI)

#DEM
DEMpath <- "C:/Users/FGreifeneder/Documents/tmp_proc/RS2_argentina/DEM/SRTM1sec_argentina_stack.tif"
slp <- raster(DEMpath, band=1)
asp <- raster(DEMpath, band=2)
hgt <- raster(DEMpath, band=3)
slpres <- resample(slp, rs2hh)
aspres <- resample(asp, rs2hh)
hgtres <- resample(hgt, rs2hh)
remove(slp, asp, hgt)

#------------------------------------------------------------
#Estimate soil moisture
#------------------------------------------------------------

smc_map <- raster(rs2hh)
# load SVR model
load("./Radarsat2_muliple_testsites/ML_model_all.dat")

#pb <- txtProgressBar(min = 0, max = nrow(rs2hh)*ncol(rs2hh), style=2)
#cntr <- 1

cl <- makeCluster(4)
registerDoSNOW(cl)

data <- foreach(i = 1:(nrow(rs2hh)*ncol(rs2hh)), .combine = cbind) %dopar% {

    if (rs2hh[i] != 0){
      # extract data
      tmpdf <- data.frame(hh = 10*log10(rs2hh[i]),
                          hv = 10*log10(rs2hv[i]),
                          lia = rs2lia[i],
                          ndvi = NDVIres[i],
                          hgt = hgtres[i],
                          slp = slpres[i],
                          asp = aspres[i])
      # estimate smc
      smc <- predict(mod, tmpdf)
    } else {
      smc <- -1
    }
    # return
    smc
}

stopCluster(cl)
#close(pb)

smc_matrix <- matrix(data=data, nrow=nrow(rs2hh), ncol=ncol(rs2hh), byrow=TRUE)
smc_map[,] <- smc_matrix

writeRaster(smc_map, 'C:/Users/FGreifeneder/Documents/tmp_proc/RS2_argentina/processing/smc_map20150809_MLall.tif', 'GTiff')
file.remove(rs2hhpath)
file.remove(rs2hvpath)
file.remove(rs2liapath)


