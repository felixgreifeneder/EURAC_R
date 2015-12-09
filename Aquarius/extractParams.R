#extract SVR training parameters from the AQUARIUS global data set and save to file

source("Aquarius/readAquInfo.R")
source("Aquarius/write_data_to_textfile.R")
library("scales")
library("ggplot2")
library("e1071")
library("caret")

extractParams <- function(AquaPath, outpath){
  
  #get the file list
  files <- list.files(AquaPath, pattern="[[:digit:]].nc", full.names=T)
  
  #initialise data frame for results. Maximum extent 360*180 deg = 64800
  dataCollection <- data.frame(time = numeric(),
                               lat = numeric(),
                               lon = numeric(),
                               ascatF = numeric(),
                               ascatincF = numeric(),
                               ascatM = numeric(),
                               ascatincM = numeric(),
                               ascatA = numeric(),
                               ascatincA = numeric(),
                               slope = numeric(),
                               aquVV = numeric(),
                               aquVH = numeric(),
                               aquHH = numeric(),
                               aquINC = numeric(),
                               aquSFLAG = numeric(),
                               eraSMC = numeric(),
                               TEMP = numeric())
  stCounter <- 0
  
  
  
  print("Extracting data from NetCDF files ...")
  for (filepath in files){
    #iterate through each file in the list. and extract data fulfilling the requirements below.
    
    aqu_params <- readAqu(filepath)
    if (length(aqu_params) == 0) next
    stations <- unique(aqu_params$locId)
    
    #iterate through all stations
    for (stId in stations){
      stCounter <- stCounter + 1
      stIndices <- which(aqu_params$locId == stId)
      
      #extract data for the current point
      stLat <- aqu_params$lat[stId+1]
      stLon <- aqu_params$lon[stId+1]
      
      time <- aqu_params$time[stIndices]
      
      ascatF <- aqu_params$ascat_sigf[stIndices]
      ascatincF <- aqu_params$ascat_incf[stIndices]
      ascatM <- aqu_params$ascat_sigm[stIndices]
      ascatincM <- aqu_params$ascat_incm[stIndices]
      ascatA <- aqu_params$ascat_siga[stIndices]
      ascatincA <- aqu_params$ascat_inca[stIndices]
      slope <- aqu_params$slope[stIndices]
      
      aquVV <- aqu_params$aqu_vv[stIndices]
      aquVH <- aqu_params$aqu_vh[stIndices]
      aquHH <- aqu_params$aqu_hh[stIndices]
      aquINC <- aqu_params$aqu_inc[stIndices]
      aquSFLAG <- aqu_params$aqu_scat_flag[stIndices]
      
      eraSMC <- aqu_params$smc[stIndices]
      TEMP <- aqu_params$soiltemp[stIndices]
      ICE <- aqu_params$ice[stIndices]
      eraSNOW <- aqu_params$snow[stIndices]
      vhvv <- 10^(aquVH/10)/10^(aquVV/10) 
      
      val <- which(eraSMC <= 0.5 & eraSMC >= 0 &  
                     TEMP > 4 & 
                     ICE == 0 &
                     eraSNOW == 0 &
                     vhvv > 0.1 &
                     vhvv < 0.4 &
                     aquSFLAG < 268435456)
      #print(paste(length(val),"/",length(aquSFLAG),sep=""))
#       if (aquSFLAG == 2147483648 |
#           aquSFLAG == 1073741824 | 
#           aquSFLAG == 536870912 |
#           aquSFLAG == 268435456 |
#           aquSFLAG == 2415919104 |
#           aquSFLAG == 2684354560 |
#           aquSFLAG == 3221225472 |
#           aquSFLAG == 1610612736 |
#           aquSFLAG == 1342177280 |
#           aquSFLAG == 805306368 |
#           aquSFLAG == 3758096384 |
#           aquSFLAG == 3489660928 |
#           aquSFLAG == 2952790016 |
#           aquSFLAG == 1879048192 |
#           aquSFLAG == 4026531840) print("RFI")
      
      if (length(val) > 0) {
        #assemble data in data frame 
        dataCollection <- rbind(dataCollection, cbind(rep(stLat, length(val)), rep(stLon, length(val)), 
                                                      time[val],
                                                      ascatF[val],
                                                      ascatincF[val],
                                                      ascatM[val],
                                                      ascatincM[val],
                                                      ascatA[val],
                                                      ascatincA[val],
                                                      slope[val],
                                                      aquVV[val],
                                                      aquVH[val],
                                                      aquHH[val],
                                                      aquINC[val],
                                                      aquSFLAG[val],
                                                      eraSMC[val],
                                                      TEMP[val]))
        
      } else {
        stCounter <- stCounter-1
      }
      
    }
    
  }
  
  names(dataCollection) <- c("lat", "lon", "time", "ascatF", "ascatincF", "ascatM", "ascatincM", "ascatA", "ascatincA", "slope",
                             "aquVV", "aquVH", "aquHH", "aquINC", "aquSFLAG", "smc", "TEMP")
  
  #plot global grid with valid points
  map <- borders("world", colour="gray50", fill="gray50")
  dev.new(width=10, height=6, noRStudioGD=T)
  gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat), data=dataCollection, colour="red", size=1.2)
  gridPlt + ggtitle(paste("Global grid - Beam", substr(outpath, 52, 52), sep=""))
  ggsave(file=paste(outpath, "global_grid_flt.png", sep=""), dpi=300)
  #graphics.off()
  
  
  #save data to binary file
  save(dataCollection, file=paste(outpath, "parameters_global_R", sep=""))
  writetxtData(dataCollection, paste(outpath, "parameters_global.txt", sep=""))
  
}