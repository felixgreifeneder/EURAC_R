#this routine computes correlations between smc and vv or vh, respectively

source("Aquarius/readAquInfo.R")
library("scales")
library("ggplot2")
library("RColorBrewer")
library("hexbin")

#get the file list
files <- list.files("X:/Workspaces/GrF/01_Data/Aquarius/aquarius_global/beam1/", pattern="[[:digit:]].nc", full.names=T)

#initialise data frame for results. Maximum extent 360*180 deg = 64800
SMCdata <- data.frame(stVV = numeric(), stVH = numeric(), stSMC = numeric(), stERASMC <- numeric())
stCounter <- 0

for (filepath in files){
  #iterate through each file in the list. Calculate correlation between smc and vv/vh
  
  aqu_params <- readAqu(filepath)
  if (length(aqu_params) == 0) next
  stations <- unique(aqu_params$locId)
  
  #iterate through all stations
  for (stId in stations){
    stCounter <- stCounter + 1
    stIndices <- which(aqu_params$locId == stId)
    
    stLat <- aqu_params$lat[stId+1]
    stLon <- aqu_params$lon[stId+1]
    stVV <- aqu_params$vv[stIndices]
    stVH <- aqu_params$vh[stIndices]
    stSMC <- aqu_params$smc[stIndices]
    stTMP <- aqu_params$surftmp[stIndices]
    stERASMC <- aqu_params$erasmc[stIndices]
    stERATMP <- aqu_params$surftmpera[stIndices]
    
    stIceflag <- aqu_params$iceflag[stIndices]
    
    val <- which(stIceflag == 0 & stSMC >= 0 & stSMC <= 1 & stTMP > 273.15 & stERATMP > 0)
    if (length(val) == 0) {
      stCounter <- stCounter-1
      next
    }
    
    data <- cbind(stVV[val], stVH[val], stSMC[val], stERASMC[val])
    SMCdata <- rbind(SMCdata, data)
    
    rm(stLat,stLon,stVV,stVH,stSMC,stERASMC,stIceflag)
    
  }
  
}

rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)

dev.new(width=8, height=6, noRStudioGD=T)
hexbinplot(V1~V3, data=SMCdata, colramp=rf, aspect=1, type="r", trans=log, inv=exp, main="SMC vs VV")

dev.new(width=8, height=6, noRStudioGD=T)
hexbinplot(V2~V3, data=SMCdata, colramp=rf, aspect=1, type="r", trans=log, inv=exp, main="SMC vs VH")

dev.new(width=8, height=6, noRStudioGD=T)
hexbinplot(V1~V4, data=SMCdata, colramp=rf, aspect=1, type="r", trans=log, inv=exp, main="SMC-ERA vs VV")

dev.new(width=8, height=6, noRStudioGD=T)
hexbinplot(V2~V4, data=SMCdata, colramp=rf, aspect=1, type="r", trans=log, inv=exp, main="SMC-ERA vs VH")

