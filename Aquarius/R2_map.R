#this routine computes correlations between smc and vv or vh, respectively

source("Aquarius/readAquInfo.R")
library("scales")
library("ggplot2")

#get the file list
files <- list.files("X:/ProjectData/EUMETSAT_SCA Cross-pol/TUW_aquarius_global_v2/beam2/", pattern="[[:digit:]].nc", full.names=T)

#initialise data frame for results. Maximum extent 360*180 deg = 64800
R2data <- data.frame(x = rep(-9999, times=64800), 
                     y = rep(-9999, times=64800), 
                     R2ascatF = rep(-9999, times=64800),
                     R2ascatM = rep(-9999, times=64800),
                     R2ascatA = rep(-9999, times=64800),
                     R2aquvv = rep(-9999, times=64800), 
                     R2aquvh = rep(-9999, times=64800))
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
    
    ascatF <- aqu_params$ascat_sigf[stIndices]
    ascatM <- aqu_params$ascat_sigm[stIndices]
    ascatA <- aqu_params$ascat_siga[stIndices]
    
    aquVV <- aqu_params$aqu_vv[stIndices]
    aquVH <- aqu_params$aqu_vh[stIndices]
    aquRFIV <- aqu_params$aqu_rfi_v[stIndices]
    aquRFIH <- aqu_params$aqu_rfi_h[stIndices]
    
    eraSMC <- aqu_params$smc[stIndices]
    TEMP <- aqu_params$soiltemp[stIndices]
    
    valASCAT <- which(eraSMC >= 0 & eraSMC <= 1 & TEMP > 0)
    
    if (length(valASCAT) > 0) {
      Corr_ascatF <- cor(eraSMC[valASCAT], y=ascatF[valASCAT], method="pearson")
      Corr_ascatM <- cor(eraSMC[valASCAT], y=ascatM[valASCAT], method="pearson")
      Corr_ascatA <- cor(eraSMC[valASCAT], y=ascatA[valASCAT], method="pearson")
    } else {
      Corr_ascatF <- -9999
      Corr_ascatM <- -9999
      Corr_ascatA <- -9999
    }
    
    valAQUVV <- which(eraSMC >= 0 & eraSMC <= 1 & TEMP > 0 & aquRFIV != 0)
    valAQUVH <- which(eraSMC >= 0 & eraSMC <= 1 & TEMP > 0 & aquRFIV != 0 & aquRFIH != 0)
    
    if (length(valAQUVV) > 0) {
      Corr_aquVV <- cor(eraSMC[valAQUVV], y=aquVV[valAQUVV], method="pearson")
    } else {
      Corr_aquVV <- -9999
    }
    
    if (length(valAQUVH) > 0) {
      Corr_aquVH <- cor(eraSMC[valAQUVH], y=aquVH[valAQUVH], method="pearson")
    } else {
      Corr_aquVH <- -9999
    }
    
    if (length(valASCAT) == 0 & length(valAQUVV) == 0 & length(valAQUVH) == 0){
      stCounter <- stCounter-1
      next
    }
    R2data$x[stCounter] <- stLon
    R2data$y[stCounter] <- stLat
    R2data$R2ascatF[stCounter] <- Corr_ascatF
    R2data$R2ascatM[stCounter] <- Corr_ascatM
    R2data$R2ascatA[stCounter] <- Corr_ascatA
    R2data$R2aquvv[stCounter] <- Corr_aquVV
    R2data$R2aquvh[stCounter] <- Corr_aquVH
    
    rm(stLat, stLon, ascatF, ascatM, ascatA, aquVV, aquVH, eraSMC, TEMP)
    
  }
  
}

map <- borders("world", colour="gray50", fill="gray50")

val <- which(R2data$R2ascatF != -9999)
dev.new(width=10, height=6, noRStudioGD=T)
ascatFplt <- ggplot() + map + geom_point(aes(x=x, y=y, colour=R2ascatF), data=R2data[val,], size=1.2) + scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"))
ascatFplt + ggtitle("Correlation SMC vs ASCAT fore beam")

val <- which(R2data$R2ascatM != -9999)
dev.new(width=10, height=6, noRStudioGD=T)
ascatFplt <- ggplot() + map + geom_point(aes(x=x, y=y, colour=R2ascatM), data=R2data[val,], size=1.2) + scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"))
ascatFplt + ggtitle("Correlation SMC vs ASCAT mid beam")

val <- which(R2data$R2ascatA != -9999)
dev.new(width=10, height=6, noRStudioGD=T)
ascatFplt <- ggplot() + map + geom_point(aes(x=x, y=y, colour=R2ascatA), data=R2data[val,], size=1.2) + scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"))
ascatFplt + ggtitle("Correlation SMC vs ASCAT aft beam")

val <- which(R2data$R2aquvv != -9999)
dev.new(width=10, height=6, noRStudioGD=T)
ascatFplt <- ggplot() + map + geom_point(aes(x=x, y=y, colour=R2aquvv), data=R2data[val,], size=1.2) + scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"))
ascatFplt + ggtitle("Correlation SMC vs Aquarius VV")

val <- which(R2data$R2aquvh != -9999)
dev.new(width=10, height=6, noRStudioGD=T)
ascatFplt <- ggplot() + map + geom_point(aes(x=x, y=y, colour=R2aquvh), data=R2data[val,], size=1.2) + scale_colour_gradient2(low=muted("red"), mid="white", high=muted("blue"))
ascatFplt + ggtitle("Correlation SMC vs Aquarius VH")

