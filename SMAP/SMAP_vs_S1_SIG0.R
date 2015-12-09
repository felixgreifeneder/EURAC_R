#plot the correlation between one scene of SMAP SIG0 vs on scene of S1 SIG0 (only grasslands)

library(raster)
library(ggplot2)
library(e1071)

S1vv <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAPvsS1/S1A_IW_GRDH_1SDV_20150626T052642_20150626T052711_006540_008B28_B9D9_VV.tif")
S1vh <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAPvsS1/S1A_IW_GRDH_1SDV_20150626T052642_20150626T052711_006540_008B28_B9D9_VH.tif")
SMAPvv_fore <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/L1C_SIG0/tifs/SMAP_L1C_S0_HIRES_02145_D_20150627T054956_R11850_001_dB_vvfore.tif")
SMAPvv_aft <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/L1C_SIG0/tifs/SMAP_L1C_S0_HIRES_02145_D_20150627T054956_R11850_001_dB_vvaft.tif")
SMAPxpol_fore <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/L1C_SIG0/tifs/SMAP_L1C_S0_HIRES_02145_D_20150627T054956_R11850_001_dB_vhfore.tif")
SMAPxpol_aft <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAP/L1C_SIG0/tifs/SMAP_L1C_S0_HIRES_02145_D_20150627T054956_R11850_001_dB_vhaft.tif")

S1vv[S1vv > -2] <- NA
S1vh[S1vh > -2] <- NA

#S1vv <- projectRaster(S1vv, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", method="ngb")
#S1vh <- projectRaster(S1vh, crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs", method="ngb")

SMAPvv_fore <- crop(SMAPvv_fore, S1vv)
SMAPvv_aft <- crop(SMAPvv_aft, S1vv)
SMAPxpol_fore <- crop(SMAPxpol_fore, S1vh)
SMAPxpol_aft <- crop(SMAPxpol_aft, S1vh)

S1vv_res <- 10*log10(resample(10^(S1vv/10), SMAPvv_fore, method="bilinear"))
S1vh_res <- 10*log10(resample(10^(S1vh/10), SMAPxpol_fore, method="bilinear"))

SMAPvv_fore <- mask(SMAPvv_fore, S1vv_res)
SMAPvv_aft <- mask(SMAPvv_aft, S1vv_res)
SMAPxpol_fore <- mask(SMAPxpol_fore, S1vh_res)
SMAPxpol_aft <- mask(SMAPxpol_aft, S1vh_res)

GC <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/SMAPvsS1/GC_res.tif")
GC <- resample(GC, SMAPvv_fore, method="ngb")
#GC[GC!=14 & GC!=20 & GC!=30] <- NA
GC[GC != 14 & GC != 30] <- NA

SMAPvv_fore <- mask(SMAPvv_fore, GC)
SMAPvv_aft <- mask(SMAPvv_aft, GC)
SMAPxpol_fore <- mask(SMAPxpol_fore, GC)
SMAPxpol_aft <- mask(SMAPxpol_aft, GC)
S1vv_res <- mask(S1vv_res, GC)
S1vh_res <- mask(S1vh_res, GC)

dat <- data.frame(SMAP=c((getValues(SMAPvv_fore)+getValues(SMAPvv_aft))/2,(getValues(SMAPxpol_fore)+getValues(SMAPxpol_aft))/2),
                  S1=c(getValues(S1vv_res),getValues(S1vh_res)),
                  pol=c(rep("VV", length(getValues(S1vv_res))), rep("VH", length(getValues(S1vh_res)))))


ggplot(dat, aes(x=SMAP, y=S1, colour=pol)) + 
  geom_point(shape=1, size=2) +
  theme(aspect.ratio=1) +
  scale_x_continuous(name="\nSIG0 SMAP [dB]", limits=c(-25,0)) +
  scale_y_continuous(name="SIG0 S1 [dB]\n", limits=c(-25,0)) + 
  scale_colour_hue(l=50, guide=F) +
  facet_grid(.~ pol) +
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14, face="bold")) +
  geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1)


#regression with SVR ------------------------------------------------------------------------------

dat2 <- data.frame(S1vv=getValues(S1vv_res), 
                   SMAPvv_f=getValues(SMAPvv_fore),
                   SMAPvv_a=getValues(SMAPvv_aft),
                   SMAPxpol_f=getValues(SMAPxpol_fore),
                   SMAPxpol_a=getValues(SMAPxpol_aft))

dat2 <- dat2[which(is.finite(dat2$S1vv) & is.finite(dat2$SMAPvv_f) & is.finite(dat2$SMAPvv_a) & is.finite(dat2$SMAPxpol_f) & is.finite(dat2$SMAPxpol_a)),] 

#--------------------------------------------------------------------------------------------------
#Building the SVR model
#package e1071

#setting the tuning options
SVRtuningsettings <- tune.control(sampling = "cross", 
                                  sampling.aggregate = mean,
                                  sampling.dispersion = sd, 
                                  cross = 5,
                                  best.model = TRUE,
                                  performances = TRUE)

#tuning the model (finding the best set  hyper parameters)
SVRtuning <- tune(svm, S1vv ~ SMAPvv_f + SMAPvv_a + SMAPxpol_f + SMAPxpol_a, 
                  data=dat2, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel <- SVRtuning$best.model

SMCpredicted <- predict(tunedModel,dat2)
error <- sqrt(mean((dat2$S1vv - SMCpredicted)^2))
r2 <- cor(dat2$S1vv, y=SMCpredicted, method="pearson")^2
print(paste("Error:",error,"R2:",r2))
