library(ncdf4)
library(R.matlab)

acqu <- list()
acqu[[1]] <- nc_open("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/roi_2/aquarius_beam_1.nc")
acqu[[2]] <- nc_open("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/roi_2/aquarius_beam_2.nc")
acqu[[3]] <- nc_open("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/roi_2/aquarius_beam_3.nc")

modis <- nc_open("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/roi_2/modis_vi.nc")
om_list <- list()

for (i in c(1:3)){
  print(paste("Beam", as.character(i), sep=" "))
  #find closest modis point
  acqu_lat <- ncvar_get(acqu[[i]], varid="lat")
  acqu_lon <- ncvar_get(acqu[[i]], varid="lon")
  modis_lat <- ncvar_get(modis, varid="lat")
  modis_lon <- ncvar_get(modis, varid="lon")
  
#   dist <- rep(0, length(modis_lat))
#   for (loc_ind in c(1:length(modis_lat))){
#     dist[loc_ind] <- sqrt((modis_lat[loc_ind]-acqu_lat)^2+(modis_lon[loc_ind]-acqu_lon)^2)
#   }
#   min_dist_loc <- which.min(dist)
  
  #retrieve variables
  acqu_time <- ncvar_get(acqu[[i]], varid="time")
  acqu_sm <- ncvar_get(acqu[[i]], varid="anc_sm")
  acqu_vv <- ncvar_get(acqu[[i]], varid="vv_toa")
  acqu_vh <- ncvar_get(acqu[[i]], varid="vh_toa")
  acqu_inc <- ncvar_get(acqu[[i]], varid="inc")
  modis_time <- ncvar_get(modis, varid="time")
  modis_ndvi <- ncvar_get(modis, varid="ndvi")
  
  out_matrix <- matrix(0, length(acqu_sm), 6)
  out_matrix[,1] <- acqu_time
  out_matrix[,2] <- acqu_sm
  out_matrix[,3] <- acqu_vv
  out_matrix[,4] <- acqu_vh
  out_matrix[,5] <- acqu_inc
  
  #for each acquarius time point, find the nearast modis acquisition
  for (ac_ind in c(1:length(acqu_sm))){
    dT <- rep(0, nrow(modis_ndvi))
    for (mod_ind in c(1:length(modis_ndvi))){
      dT[mod_ind] <- abs(acqu_time[ac_ind]-modis_time[mod_ind])
    }
    min_dT <- which.min(dT)
    out_matrix[ac_ind,6] <- mean(modis_ndvi[min_dT, ])
  }
  
  om_list[[i]] <- out_matrix
}

matVect <- rbind(om_list[[1]], om_list[[2]], om_list[[3]])
writeMat("X:/Workspaces/GrF/Processing/EUMETSAT/SVR/ROI2/params_roi2.mat", matVect=matVect)

