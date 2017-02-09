library(raster)
library(ggplot2)
library(maptools)

landsat <- stack("//projectdata.eurac.edu//projects//ActivityData//Padovano_PhD//PaA//08_SMEX02//landsat//20020708_L7//20020708_L7//01_stk//20020708_L7_stk")
NAvalue(landsat) <- -9999

c.hh <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=1)
c.vv <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=2)
c.hv <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=3)
c.vh <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=4)
l.hh <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=5)
l.vv <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=6)
l.hv <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=7)
l.vh <- raster("//projectdata.eurac.edu/projects/ActivityData/Padovano_PhD/PaA/08_SMEX02/AiRSAR data decomp/20020708_sar/03_geocodificato/8_July_brk_geo.tif", band=8)

NAvalue(c.hh) <- 0
NAvalue(c.vv) <- 0
NAvalue(c.hv) <- 0
NAvalue(c.vh) <- 0
NAvalue(l.hh) <- 0
NAvalue(l.vv) <- 0
NAvalue(l.hv) <- 0
NAvalue(l.vh) <- 0

# resample
landsat <- resample(landsat, c.vv)


# masking
mask.shape <- readShapePoly("C:/Users/FGreifeneder/Documents/tmp_proc/ms0/AIRSAR/maskUTM2.shp", 
                            proj4string = CRS("+proj=utm +zone=15 +datum=WGS84 +units=m +no_defs"))

# compute polarization ratios and ndvi
#c.vvvh <- c.vv/c.vh
#l.vvvh <- l.vv/l.vh
c.vvvh <- (8*c.hv)/(c.hh+c.vv+(2*c.hv))
l.vvvh <- (8*l.hv)/(l.hh+l.vv+(2*l.hv))
ndvi <- (raster(landsat, layer=4)-raster(landsat, layer=3)) / (raster(landsat, layer=4) + raster(landsat, layer=3))

# create layer stack
datastack <- stack(c.vv,c.hh,c.hv,c.vh,l.vv,l.hh,l.hv,l.vh, c.vvvh, l.vvvh, ndvi)

# extract samples
extracted.samples <- extract(datastack, mask.shape, fun=mean, na.rm=T, layer=1, nl=11)
colnames(extracted.samples) <- c("c.vv","c.hh","c.hv","c.vh","l.vv","l.hh","l.hv","l.vh","c.vvvh","l.vvvh","ndvi")
extracted.samples <- data.frame(extracted.samples)

# mask values
#extracted.samples <- extracted.samples[which(extracted.samples$c.vvvh > 1),]

# compute fit 
linear.fit <- lm(c.vvvh ~ l.vvvh, data=extracted.samples)
linear.coefficients <- coefficients(linear.fit)
c.vvvh.pred <- linear.coefficients[1] + linear.coefficients[2]*extracted.samples$l.vvvh
# correlation
cor(extracted.samples$c.vvvh, c.vvvh.pred)
# rmse
sqrt(sum((c.vvvh.pred-extracted.samples$c.vvvh)^2))
predicted.samples <- data.frame(c.vvvh=extracted.samples$c.vvvh, c.vvvh.pred=c.vvvh.pred)

# c.vvvh vs l.vvvh
ggplot(extracted.samples, aes(x=10*log10(l.vvvh), y=10*log10(c.vvvh))) + 
  geom_point(shape=1) + 
  geom_smooth(method='lm') + 
  xlim(-3,5) + 
  ylim(-3,5) +
  xlab('L-Band VV/VH [dB]') +
  ylab('C-Band VV/VH [dB]') +
  theme(axis.title=element_text(size=10)) + 
  annotate("text", x=1, y=-1, label="y = 0.24x + 1.33")

# predicted c.vvvh vs c.vvvh
ggplot(predicted.samples, aes(x=10*log10(c.vvvh.pred), y=10*log10(c.vvvh))) + 
  geom_point(shape=1) + 
  geom_smooth(method='lm') + 
  xlim(0,4) + 
  ylim(0,4) + 
  xlab('Predicted C-Band VV/VH [dB]') + 
  ylab('C-Band VV/VH') +
  theme(axis.title=element_text(size=10)) +
  annotate("text", x=1, y=0.5, label="R: 0.75\nRMSE: 0.79")

ggplot(extracted.samples, aes(x=ndvi, y=10*log10(c.vvvh))) + geom_point(shape=1)
ggplot(extracted.samples, aes(x=ndvi, y=10*log10(l.vvvh))) + geom_point(shape=1)

