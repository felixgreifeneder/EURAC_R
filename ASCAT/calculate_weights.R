#calculates the station scaling weights for the ascat footprint

library("raster")

calcWeights <- function(path){
  
  class_raster <- raster(path)
  class_image <- getValues(class_raster)
  
  nb1 <- length(which(class_image == 6))
  nb2 <- length(which(class_image == 5))
  nm1 <- length(which(class_image == 4))
  nm2 <- length(which(class_image == 8))
  nm3 <- length(which(class_image == 7))
  nm4 <- length(which(class_image == 3))
  nm5 <- length(which(class_image == 2))
  ns4 <- length(which(class_image == 1))
  
  nrPixels <- nb1+nb2+nm1+nm2+nm3+nm4+nm5+ns4
  
  wb1 <- nb1/nrPixels
  wb2 <- nb2/nrPixels
  wm1 <- nm1/nrPixels
  wm2 <- nm2/nrPixels
  wm3 <- nm3/nrPixels
  wm4 <- nm4/nrPixels
  wm5 <- nm5/nrPixels
  ws4 <- ns4/nrPixels
  
  print(paste("Wb1:", wb1))
  print(paste("Wb2:", wb2))
  print(paste("Wm1:", wm1))
  print(paste("Wm2:", wm2))
  print(paste("Wm3:", wm3))
  print(paste("Wm4:", wm4))
  print(paste("Wm5:", wm5))
  print(paste("Ws4:", ws4))
  
  weights <- data.frame(WB1=wb1,WB2=wb2,WM1=wm1,WM2=wm2,WM3=wm3,WM4=wm4,WM5=wm5,WS4=ws4)
  
  return(weights)
  
  
}