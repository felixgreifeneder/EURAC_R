#run exec.R first

library("ggplot2")

nRows <- length(which(B1SIG0$vv_fore != -9999 & !(B1SIG0$rflag >= 8 & B1SIG0$rflag < 16) & B1SIG0$dist < 0.0045))*6 +
         length(which(B2SIG0$vv_fore != -9999 & !(B2SIG0$rflag >= 8 & B2SIG0$rflag < 16) & B2SIG0$dist < 0.0045))*6 +
         length(which(B3SIG0$vv_fore != -9999 & !(B3SIG0$rflag >= 8 & B3SIG0$rflag < 16) & B3SIG0$dist < 0.0045))*6

Stnames <- names(B1SIG0)

Df <- data.frame(SIG0=numeric(length = nRows), 
                  elev_std=numeric(length = nRows),
                  pol=character(length=nRows),
                  station=character(length=nRows), stringsAsFactors = F)


#fill all B1 data into the data frame -------------------------------------------------------------
cntr <- 1
vRows <- which(B1SIG0$vv_fore != -9999 & !(B1SIG0$rflag >= 8 & B1SIG0$rflag < 16) & B1SIG0$dist < 0.0045)

for (i in 1:6){
  
  for (j in vRows){
    
    Df$SIG0[cntr] <- B1SIG0[j,i]
    Df$elev_std[cntr] <- B1SIG0$elev_std[j]
    Df$pol[cntr] <- Stnames[i]
    Df$station[cntr] <- "B1000"
    cntr <- cntr + 1
    
  }
  
}

#B2------------------------------------------------------------------------------------------------
vRows <- which(B2SIG0$vv_fore != -9999 & !(B2SIG0$rflag >= 8 & B2SIG0$rflag < 16) & B2SIG0$dist < 0.0045)

for (i in 1:6){
  
  for (j in vRows){
    
    Df$SIG0[cntr] <- B2SIG0[j,i]
    Df$elev_std[cntr] <- B2SIG0$elev_std[j]
    Df$pol[cntr] <- Stnames[i]
    Df$station[cntr] <- "B1500"
    cntr <- cntr + 1
    
  }
  
}

#B3------------------------------------------------------------------------------------------------
vRows <- which(B3SIG0$vv_fore != -9999 & !(B3SIG0$rflag >= 8 & B3SIG0$rflag < 16) & B3SIG0$dist < 0.0045)

for (i in 1:6){
  
  for (j in vRows){
    
    Df$SIG0[cntr] <- B3SIG0[j,i]
    Df$elev_std[cntr] <- B3SIG0$elev_std[j]
    Df$pol[cntr] <- Stnames[i]
    Df$station[cntr] <- "B2000"
    cntr <- cntr + 1
    
  }
  
}

#Plot----------------------------------------------------------------------------------------------

dev.new()

ggplot(data=Df, aes(x=elev_std,y=10*log10(SIG0))) + 
  geom_point()  +
  theme(aspect.ratio=1) +
  geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1) +
  facet_grid(pol ~ station) + ggtitle("Elev. Std-Dev vs SIG0 \n")



