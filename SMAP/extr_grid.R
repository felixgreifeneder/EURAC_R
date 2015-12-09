#this routine extracts grid-points within the area of interest for grid definition
library(rhdf5)


h5list <- list.files("D:/SMAP/L1C_SIG0", pattern=".h5$", full.names=TRUE)

xylist <- data.frame(lat=as.numeric(), lon=as.numeric())
valid_files <- 0
posCntr <- 1

total <- length(h5list)
pb <- txtProgressBar(min=0, max=total, style=3)

for (find in 1:length(h5list)){
  
  latlist <- h5read(h5list[find], "/Sigma0_Data/cell_lat")
  lonlist <- h5read(h5list[find], "/Sigma0_Data/cell_lon")
  
  aoiIdx <- which(latlist >= 44 & latlist <= 50 &
                  lonlist >= 4.5 & lonlist <= 17.5)
  
  if (length(aoiIdx)== 0) {
    H5close()
    setTxtProgressBar(pb, find)
    next
  }
  # if (valid_files == 0) {
    
    #for (i in 1:length(aoiIdx)){
      
#       xylist[posCntr,1] <- latlist[aoiIdx[i]]
#       xylist[posCntr,2] <- lonlist[aoiIdx[i]]
      xylist <- rbind(xylist, data.frame(lat=latlist[aoiIdx], lon=lonlist[aoiIdx]))
      posCntr <- posCntr+1
      
    #} 
    
    valid_files <- valid_files+1
    
#   } else {
#     
#     for (i in 1:length(aoiIdx)){
#       
#       #if (any(xylist[,1]==latlist[aoiIdx[i]] & xylist[,2]==lonlist[aoiIdx[i]]) == F){
#         xylist[posCntr,1] <- latlist[aoiIdx[i]]
#         xylist[posCntr,2] <- lonlist[aoiIdx[i]]
#         posCntr <- posCntr+1
#       #}
#       
#     }
#     
#   }
  
  H5close()
  setTxtProgressBar(pb, find)
  
}