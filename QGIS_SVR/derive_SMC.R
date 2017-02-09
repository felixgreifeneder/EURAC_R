#function definition for parallel SMC predictions
library(raster)
library(gbm)
library(caret)
library(ggplot2)
library(foreach)
library(doSNOW)

parallel_predictions<-function(fit,testing)
{
  cl<-makeCluster(4)
  registerDoSNOW(cl)
  num_splits<-4
  split_testing<-sort(rank(1:nrow(testing))%%4)
  predictions<-foreach(i=unique(split_testing),
                       .combine=c,.packages=c("caret")) %dopar% {
                         as.numeric(predict(fit,newdata=testing[split_testing==i,]))
                       }
  stopCluster(cl)
  return(predictions)
}

#script for the production of smc maps


#working dir
Working_Directory <- 'C:/Users/FGreifeneder/Documents/work'
SAR_SIG0HH <- "X:/Workspaces/GrF/Processing/SVR_QGIS/RS2_18112014_HH_SIG0.tif"
SAR_SIG0HV <- "X:/Workspaces/GrF/Processing/SVR_QGIS/RS2_18112014_HV_SIG0.tif"
SAR_SIG0VV <- ''
SAR_SIG0VH <- ''
LIA <- "X:/Workspaces/GrF/Processing/SVR_QGIS/RS2_18112014_LIA.tif"
DEM <- "X:/Workspaces/GrF/Processing/SVR_QGIS/DEM_stack.tif"
NDVI <- "X:/Workspaces/GrF/Processing/SVR_QGIS/NDVI_17112014.tif"
Land_Cover <- ''
Model_File <- "X:/Workspaces/GrF/Processing/SVR_QGIS/qgismodel.dat"
Output_File <- "X:/Workspaces/GrF/Processing/SVR_QGIS/smc_map.tif"


# -------------------------------------------------------------
# Initiate required datasets
# -------------------------------------------------------------

# initiate data stack
data_stack <- stack()

# SAR dataset


# copy fles to working directory, for faster data access
if (SAR_SIG0VV != '') {
  file.copy(SAR_SIG0VV, paste(Working_Directory, basename(SAR_SIG0VV), sep="/"))
  vvpath <- paste(Working_Directory, basename(SAR_SIG0VV), sep="/")
  vv <- raster(vvpath)
  data_stack <- addLayer(data_stack, vv)
}
if (SAR_SIG0VH != '') {
  file.copy(SAR_SIG0VH, paste(Working_Directory, basename(SAR_SIG0VH), sep="/"))
  vhpath <- paste(Working_Directory, basename(SAR_SIG0VH), sep="/")
  vh <- raster(vhpath)
  data_stack <- addLayer(data_stack, vh)
}
if (SAR_SIG0HH != '') {
  file.copy(SAR_SIG0HH, paste(Working_Directory, basename(SAR_SIG0HH), sep="/"))
  hhpath <- paste(Working_Directory, basename(SAR_SIG0HH), sep="/")
  hh <- raster(hhpath)
  data_stack <- addLayer(data_stack, hh)
}
if (SAR_SIG0HV != '') {
  file.copy(SAR_SIG0HV, paste(Working_Directory, basename(SAR_SIG0HV), sep="/"))
  hvpath <- paste(Working_Directory, basename(SAR_SIG0HV), sep="/")
  hv <- raster(hvpath)
  data_stack <- addLayer(data_stack, hv)
}
if (LIA != '') {
  file.copy(LIA, paste(Working_Directory, basename(LIA), sep="/"))
  liapath <- paste(Working_Directory, basename(LIA), sep="/")
  lia <- raster(liapath)
  liares <- resample(lia, data_stack)
  remove(lia)
  data_stack <- addLayer(data_stack, liares)
}


# DEM
if (DEM != '') {
  slp <- raster(DEM, band=1) # Slope
  asp <- raster(DEM, band=2) # Aspect
  hgt <- raster(DEM, band=3) # Height
  slpres <- resample(slp, hh)
  aspres <- resample(asp, hh)
  hgtres <- resample(hgt, hh)
  remove(slp, asp, hgt)
  data_stack <- addLayer(data_stack, hgtres)
  data_stack <- addLayer(data_stack, slpres)
  data_stack <- addLayer(data_stack, aspres)
}


# NDVI
if (NDVI != '') {
  NDVI <- raster(NDVI)
  NDVIres <- resample(NDVI, hh, method="bilinear")
  remove(NDVI)
  data_stack <- addLayer(data_stack, NDVIres)
}

#------------------------------------------------------------
#Estimate soil moisture
#------------------------------------------------------------

# Initialise SMC map (sampe extent/resolution as SAR input data)
smc_map <- raster(nrows=nrow(data_stack), 
                  ncols=ncol(data_stack), 
                  vals=-1, 
                  ext=extent(data_stack),
                  crs=crs(data_stack))
# load SVR model
load(Model_File)

# convert image stack to matrix
data_stack_mat <- as.matrix(data_stack)
smc_vector <- parallel_predictions(mod, data_stack_mat)
remove(data_stack_mat)
gc()

smc_matrix <- matrix(data=smc_vector, nrow=nrow(data_stack), ncol=ncol(data_stack), byrow=TRUE)
smc_map[,] <- smc_matrix

writeRaster(smc_map, Output_File, 'GTiff')
if (SAR_SIG0HH != "") {
  file.remove(hhpath)
}
if (SAR_SIG0VV != "") {
  file.remove(vvpath)
}
if (SAR_SIG0VH != "") {
  file.remove(vhpath)
}
if (SAR_SIG0HV != "") {
  file.remove(hvpath)
}
if (LIA != "") {
  file.remove(liapath)
}

