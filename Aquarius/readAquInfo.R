#read parameters from the TUW aquarius dataset: vv, vh, smc, slope, soiltmp, inc, lat, lon

library(ggmap, ggplot2)
library(ncdf4)
#library(rworldmap)

readAqu <- function(path){
  
  #error handling with tryCatch()
  output <- tryCatch({
    #open netcdf file
    aqu_file <- nc_open(path)
    
    #extract variables
    lat <- ncvar_get(aqu_file, varid="lat")
    lon <- ncvar_get(aqu_file, varid="lon")
    
    locId <- ncvar_get(aqu_file, varid="locationIndex")
    time <- ncvar_get(aqu_file, varid="time")
    smc <- ncvar_get(aqu_file, varid="era_39")
    
    ascat_sigf <- ncvar_get(aqu_file, varid="ascat_sigf")
    ascat_incf <- ncvar_get(aqu_file, varid="ascat_incf")
    ascat_sigm <- ncvar_get(aqu_file, varid="ascat_sigm")
    ascat_incm <- ncvar_get(aqu_file, varid="ascat_incm")
    ascat_siga <- ncvar_get(aqu_file, varid="ascat_siga")
    ascat_inca <- ncvar_get(aqu_file, varid="ascat_inca")
    
    slope <- ncvar_get(aqu_file, varid="warp_slop")
    
    aqu_vv <- ncvar_get(aqu_file, varid="aquarius_vv_toa")
    aqu_vh <- ncvar_get(aqu_file, varid="aquarius_vh_toa")
    aqu_hh <- ncvar_get(aqu_file, varid="aquarius_vh_toa")
    aqu_inc <- ncvar_get(aqu_file, varid="aquarius_inc")
    aqu_scat_flag <- ncvar_get(aqu_file, varid="aquarius_scat_flags")
    
    soiltemp <- ncvar_get(aqu_file, varid="era_139")
    snow <- ncvar_get(aqu_file, varid="era_141")
    ice <- ncvar_get(aqu_file, varid="aquarius_scat_ice_frac")
    
    nc_close(aqu_file)
    
    print(basename(path))
    
    f <- list(lat=lat, 
              lon=lon,
              locId=locId,
              time=time,
              smc=smc,
              slope=slope,
              ascat_sigf=ascat_sigf,
              ascat_incf=ascat_incf,
              ascat_sigm=ascat_sigm,
              ascat_incm=ascat_incm,
              ascat_siga=ascat_siga,
              ascat_inca=ascat_inca,
              aqu_vv=aqu_vv,
              aqu_vh=aqu_vh,
              aqu_hh=aqu_hh,
              aqu_inc=aqu_inc,
              aqu_scat_flag=aqu_scat_flag,
              soiltemp=soiltemp,
              snow=snow,
              ice=ice)
   
  }, warning = function(war) {
    #handling of warning
    print(war)
    print(paste("Error in file:", path))
    f <- list()
    return(f)
  }, error = function(err) {
    #handling of errors
    print(err)
    print(paste("Error in file:", path))
    f <- list()
    return(f)
  }, finally = {
  })
  
  
  
  return(output)
  
  
  #map <- get_map(location=centre, zoom=4)
  #map <- get_map()
  #map <- getMap(resolution="low")
  #mapPoints <- ggmap(map, extent="normal") + geom_point(aes(x=lon, y=lat), data = locations)
  #mapPoints <- geom_point(aes(x=lon, y=lat), data = locations)
  #plot(mapPoints)
  
  #return(output)
}