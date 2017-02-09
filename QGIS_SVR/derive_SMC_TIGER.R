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

# settings
# copy all datasets to a working directory. Recommend if no all datasets are stored locally
copy_to_work_dir <- T

# working dir
# must be defined if copy_to_work_dir <- T
Working_Directory <- "C:\\Users\\FGreifeneder\\Documents\\work"
SAR_SIG0VV <- "Q:/ESA_TIGER/code/data_and_scripts/Sigma0_VV_UTM37s_subset.tif"
PLIA <- "Q:/ESA_TIGER/code/data_and_scripts/plia_utm37s_subset.tif"
S2_b2 <- "Q:/ESA_TIGER/code/data_and_scripts/B02_subset.tif"
S2_b3 <- "Q:/ESA_TIGER/code/data_and_scripts/B03_subset.tif"
S2_b4 <- "Q:/ESA_TIGER/code/data_and_scripts/B04_subset.tif"
S2_b5 <- "Q:/ESA_TIGER/code/data_and_scripts/B05_subset.tif"
S2_b6 <- "Q:/ESA_TIGER/code/data_and_scripts/B06_subset.tif"
S2_b7 <- "Q:/ESA_TIGER/code/data_and_scripts/B07_subset.tif"
S2_b8 <- "Q:/ESA_TIGER/code/data_and_scripts/B08_subset.tif"
Land_Cover <- "Q:/ESA_TIGER/code/data_and_scripts/kajiado_2016_lulc_utm37s_raster_subset.tif"

Model_File <- "X:/Workspaces/GrF/Processing/ESA_TIGER/s1_nolia_model.dat"
Output_File <- "X:/Workspaces/GrF/Processing/ESA_TIGER/smc_map16082016_s1_nolia.tif"


# -------------------------------------------------------------
# Initiate required datasets
# -------------------------------------------------------------

# initiate data stack
data_stack <- stack()

# SAR dataset


# copy fles to working directory, for faster data access
if (SAR_SIG0VV != '') {
  if (copy_to_work_dir == T){
    file.copy(SAR_SIG0VV, paste(Working_Directory, basename(SAR_SIG0VV), sep="/"))
    vvpath <- paste(Working_Directory, basename(SAR_SIG0VV), sep="/")
  } else {
    vvpath <- SAR_SIG0VV
  }
  vv <- raster(vvpath)
  vv <- 10*log10(vv)
  names(vv) <- "Sig0_VV"
  data_stack <- addLayer(data_stack, vv)
}
if (PLIA != '') {
  if (copy_to_work_dir == T){
    file.copy(PLIA, paste(Working_Directory, basename(PLIA), sep="/"))
    liapath <- paste(Working_Directory, basename(PLIA), sep="/")
  } else {
    liapath <- PLIA
  }
  lia <- raster(liapath)
  liares <- resample(lia, data_stack)
  remove(lia)
  names(liares) <- "PLIA"
  data_stack <- addLayer(data_stack, liares)
}
# S2
if (S2_b2 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b2, paste(Working_Directory, basename(S2_b2), sep="/"))
    s2path1 <- paste(Working_Directory, basename(S2_b2), sep="/")
  } else {
    s2path1 <- S2_b2
  }
  s2 <- raster(s2path1)
  s2_b2res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b2res) <- "S_Band_2"
  data_stack <- addLayer(data_stack, s2_b2res)
}
if (S2_b3 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b3, paste(Working_Directory, basename(S2_b3), sep="/"))
    s2path2 <- paste(Working_Directory, basename(S2_b3), sep="/")
  } else {
    s2path2 <- S2_b3
  }
  s2 <- raster(s2path2)
  s2_b3res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b3res) <- "S_Band_3"
  data_stack <- addLayer(data_stack, s2_b3res)
}
if (S2_b4 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b4, paste(Working_Directory, basename(S2_b4), sep="/"))
    s2path3 <- paste(Working_Directory, basename(S2_b4), sep="/")
  } else {
    s2path3 <- S2_b4
  }
  s2 <- raster(s2path3)
  s2_b4res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b4res) <- "S_Band_4"
  data_stack <- addLayer(data_stack, s2_b4res)
}
if (S2_b5 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b5, paste(Working_Directory, basename(S2_b5), sep="/"))
    s2path4 <- paste(Working_Directory, basename(S2_b5), sep="/")
  } else {
    s2path4 <- S2_b5
  }
  s2 <- raster(s2path4)
  s2_b5res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b5res) <- "S_Band_5"
  data_stack <- addLayer(data_stack, s2_b5res)
}
if (S2_b6 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b6, paste(Working_Directory, basename(S2_b6), sep="/"))
    s2path5 <- paste(Working_Directory, basename(S2_b6), sep="/")
  } else {
    s2path5 <- S2_b6
  }
  s2 <- raster(s2path5)
  s2_b6res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b6res) <- "S_Band_6"
  data_stack <- addLayer(data_stack, s2_b6res)
}
if (S2_b7 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b7, paste(Working_Directory, basename(S2_b7), sep="/"))
    s2path6 <- paste(Working_Directory, basename(S2_b7), sep="/")
  } else {
    s2path6 <- S2_b7
  }
  s2 <- raster(s2path6)
  s2_b7res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b7res) <- "S_Band_7"
  data_stack <- addLayer(data_stack, s2_b7res)
}
if (S2_b8 != '') {
  if (copy_to_work_dir == T){
    file.copy(S2_b8, paste(Working_Directory, basename(S2_b8), sep="/"))
    s2path7 <- paste(Working_Directory, basename(S2_b8), sep="/")
  } else {
    s2path7 <- S2_b8
  }
  s2 <- raster(s2path7)
  s2_b8res <- resample(s2, data_stack)
  remove(s2)
  names(s2_b8res) <- "S_Band_8"
  data_stack <- addLayer(data_stack, s2_b8res)
}
if (Land_Cover != '') {
  if (copy_to_work_dir == T){
    file.copy(Land_Cover, paste(Working_Directory, basename(Land_Cover), sep="/"))
    lcpath <- paste(Working_Directory, basename(Land_Cover), sep="/")
  } else {
    lcpath <- Land_Cover
  }
  lc <- raster(lcpath)
  lcres <- resample(lc, data_stack)
  remove(lc)
  names(lcres) <- "LC"
  data_stack <- addLayer(data_stack, lcres)
}

gc()
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

if (copy_to_work_dir==T){
  if (SAR_SIG0VV != "") {
    file.remove(vvpath)
  }
  if (PLIA != "") {
    file.remove(liapath)
  }
  if (S2_b2 != "") {
    file.remove(s2path1)
  }
  if (S2_b3 != "") {
    file.remove(s2path2)
  }
  if (S2_b4 != "") {
    file.remove(s2path3)
  }
  if (S2_b5 != "") {
    file.remove(s2path4)
  }
  if (S2_b6 != "") {
    file.remove(s2path5)
  }
  if (S2_b7 != "") {
    file.remove(s2path6)
  }
  if (S2_b8 != "") {
    file.remove(s2path7)
  }
  if (Land_Cover != "") {
    file.remove(lcpath)
  }
}

