#this scrip fetches data from the mobile campaign in-situ measurements and radarsat2 acuisitions

library(raster)

#load data from mobile campaigns
load("./Radarsat2Mazia/HiResAlp_MobileCampaigns.RData")
#load table where each campaign day is linked to a radarsat acquisition
load("./Radarsat2_muliple_testsites/insitusarlin.dat")

workdir <- "C:/Users/FGreifeneder/Documents/work/"


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

CombDs <- data.frame(smc=numeric(), 
                     hh=numeric(),
                     hv=numeric(),
                     lia=numeric(),
                     height=numeric(),
                     slope=numeric(),
                     aspect=numeric(),
                     ndvi=numeric(),
                     lon=numeric(),
                     lat=numeric(),
                     date=character(), stringsAsFactors = F)
cntr <- 1

#iterate through all RS2 acquisitions
for (RowI in 1:nrow(inSituSARlink)){
  print(inSituSARlink$rs2[RowI])
  
  #check if matching field campaign exists
  if (inSituSARlink$insitu[RowI] == "NA") next
  
  #load RS2 image
  file.copy(inSituSARlink$paths[RowI], paste(workdir, basename(inSituSARlink$paths[RowI]), sep=""))
  RS2hh <- raster(paste(workdir, basename(inSituSARlink$paths[RowI]), sep=""), band=1)
  RS2hv <- raster(paste(workdir, basename(inSituSARlink$paths[RowI]), sep=""), band=2)
  RS2lia <- raster(paste(workdir, basename(inSituSARlink$paths[RowI]), sep=""), band=3)
  
  #load NDVI
  NDVIimg <- raster(inSituSARlink$ndvi[RowI])
  
  #iterate through field campaign measurements
  inSituList <- which(HiResAlp_MobileCampaigns$date == inSituSARlink$insitu[RowI])
  for (inSituind in inSituList){
    
    tmprow <- data.frame(smc=0.0, hh=0.0, hv=0.0, lia=0.0, height=0.0, slope=0.0, aspect=0.0, ndvi=0.0, lon=0.0, lat=0.0, date="", stringsAsFactors = F)
    
    
    #extract data
    lonlat <- matrix(c(HiResAlp_MobileCampaigns$longitude[inSituind], HiResAlp_MobileCampaigns$latitude[inSituind]), 1, 2)
    tmprow$smc <- HiResAlp_MobileCampaigns$SoilMoisture_mean[inSituind]
    tmprow$hh <- extract(RS2hh, 
                         lonlat,
                         layer=1, nl=1)
    tmprow$hv <- extract(RS2hv, 
                         lonlat,
                         layer=2, nl=1)
    tmprow$lia <- extract(RS2lia, 
                          lonlat,
                          layer=3, nl=1)
    
    tmprow$height <- extract(dem, lonlat)
    tmprow$slope <- extract(slope, lonlat)
    tmprow$aspect <- extract(aspect, lonlat)
    
    #tmprow$ndvi <- HiResAlp_MobileCampaigns$MODIS.NDVI.mean[inSituind]
    tmprow$ndvi <- extract(NDVIimg, lonlat)
    
    tmprow$lon <- lonlat[1]
    tmprow$lat <- lonlat[2]
    tmprow$date <- HiResAlp_MobileCampaigns$date[inSituind]
    
    CombDs[cntr,] <- tmprow
    cntr <- cntr + 1
    
  }
  
  file.remove(paste(workdir, basename(inSituSARlink$paths[RowI]), sep=""))
  
}

file.remove(dempath)
file.remove(aspectpath)
file.remove(slopepath)

save(CombDs, file="./Radarsat2_muliple_testsites/ComDs.dat")