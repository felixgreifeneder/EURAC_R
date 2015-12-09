library("ggplot2")
library("raster")

# read table with coord list

locs = read.table("C:/Users/FGreifeneder/Documents/tmp_proc/SCA_paper/insitu_coordlist.txt", sep=",", header=T)

#get landcover
globcover <- raster("C:/Users/FGreifeneder/Documents/tmp_proc/GLOBCOVER_L4_200901_200912_V2.3.tif")
lcs = extract(globcover, matrix(cbind(locs$lon, locs$lat), 34, 2))
locs = data.frame(locs, lc=as.factor(lcs))


#plot the locations
map <- borders("world", colour="gray50", fill="gray50")
#rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
#r <- rf(32)
dev.new(width=10, height=6, noRStudioGD=T)
gridPlt <- ggplot() + map + geom_point(aes(x=lon, y=lat, colour=lc), data=locs, size=3) +
           scale_colour_discrete(name="Land-Cover\n", 
                               labels=c("Cropland (50-70%)/vegetation (20-50%)",
                                        "Vegetation (50-70%)/cropland (50-70%)",
                                        "Closed brodleaved deciduous forest",
                                        "Closed needleleaved forest",
                                        "Open needleleaved forest",
                                        "Closed to open mixed forest",
                                        "Grassland (50-70%)/ forest or shrubland (20-50%)",
                                        "Closed to open shrubland",
                                        "Closed to open grassland",
                                        "Sparse vegetation",
                                        "Bare areas",
                                        "Water bodies", 
                                        "Permanent snow and ice"))
gridPlt + ggtitle("In-Situ data loactions\n")

dev.new(width=8, height=8, noRStudioGD=T)
ggplot(locs, aes(x=lc)) + geom_histogram() +
  scale_x_discrete(labels=c("Cropland (50-70%)/vegetation (20-50%)",
                            "Vegetation (50-70%)/cropland (50-70%)",
                            "Closed brodleaved deciduous forest",
                            "Closed needleleaved forest",
                            "Open needleleaved forest",
                            "Closed to open mixed forest",
                            "Grassland (50-70%)/ forest or shrubland (20-50%)",
                            "Closed to open shrubland",
                            "Closed to open grassland",
                            "Sparse vegetation",
                            "Bare areas",
                            "Water bodies", 
                            "Permanent snow and ice"), name="") +
  ylab("Count") +
  ggtitle("Distribution of land-cover classes\n") +
  coord_flip()


