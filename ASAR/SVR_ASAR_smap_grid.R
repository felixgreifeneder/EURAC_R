#Initialisation
#--------------------------------------------------------------------------------------------------
#loading required packages
library("e1071")
library("caret")
library("ggplot2")
library("reshape2")
library("hexbin")
library("RColorBrewer")
library("gridExtra")


#load the training-set
training_noScl <- read.table("X:/Workspaces/GrF/Processing/SMAP/svr/training_stations_sig0lcmask.txt",
                             header=T,
                             sep=",")

training_Scld <- read.table("X:/Workspaces/GrF/Processing/SMAP/svr/training_scaled_sig0lcmask.txt",
                            header=T,
                            sep=",")


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
SVRtuning <- tune(svm, SMC ~ SIG0 + height + aspect + slope + lc, 
                  data=training_noScl, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel_noScl <- SVRtuning$best.model

SVRtuning <- tune(svm, SMC ~ SIG0 + height + aspect + slope + lc, 
                  data=training_Scld, 
                  kernel="radial",
                  ranges = list(epsilon = 10^seq(-2,-0.5,len=3), gamma=10^seq(-2,1,len=3), cost=10^seq(-2,2,len=3)),
                  tunecontrol=SVRtuningsettings)

tunedModel_Scld <- SVRtuning$best.mode

#--------------------------------------------------------------------------------------------------
#Estimating SMC

print("Overall performance (no scaling) based on training:")
SMCpredicted_noScl <- predict(tunedModel_noScl, training_noScl)
error_noScl <- sqrt(mean((training_noScl$SMC - SMCpredicted_noScl)^2))
r2_noScl <- cor(training_noScl$SMC, y=SMCpredicted_noScl, method="pearson")^2
print(paste("Error:",error_noScl,"R2:",r2_noScl))

print("Overall performance (scaled) based on training:")
SMCpredicted_Scld <- predict(tunedModel_Scld, training_Scld)
error_Scld <- sqrt(mean((training_Scld$SMC - SMCpredicted_Scld)^2))
r2_Scld <- cor(training_Scld$SMC, y=SMCpredicted_Scld, method="pearson")^2
print(paste("Error:",error_Scld,"R2:",r2_Scld))

#Plot true vs estimated SMC
tmp <- data.frame(x=c(training_noScl$SMC, training_Scld$SMC),
                  y=c(SMCpredicted_noScl, SMCpredicted_Scld), 
                  Scaling=c(rep("no", length(training_noScl$SMC)), rep("yes", length(training_Scld$SMC))), stringsAsFactors = F)


rf <- colorRampPalette(rev(brewer.pal(11,"Spectral")))
r <- rf(32)
ggplot(tmp, aes(x=x, y=y, colour=Scaling)) + 
  geom_point(shape=1, size=3) +
  theme(aspect.ratio=1) +
  scale_x_continuous(name="\nTRUE SMC [m3m-3]", limits=c(0.1,0.5)) +
  scale_y_continuous(name="PREDICTED SMC [m3m-3]\n", limits=c(0.1,0.5)) + 
  scale_colour_hue(l=50, guide=F) +
  facet_grid(.~ Scaling) +
  theme(axis.title.x = element_text(face="bold", size=18),
        axis.text.x = element_text(size=14),
        axis.title.y = element_text(face="bold", size=18),
        axis.text.y = element_text(size=14),
        strip.text.x = element_text(size=14, face="bold")) +
  geom_smooth(method=glm, se=F, fullrange=T, linetype="dashed", size=1)
ggsave("X:/Workspaces/GrF/02_Documents/VZJ_paper/fig_for_upload/fig8_new.tiff", dpi=300)


