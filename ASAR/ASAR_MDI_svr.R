library(kernlab)
library(date)
library(myTools)
library(raster)
library(e1071)

#specify paths:
  
asar_path <- "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/ASARWS/"
model_path <- "C:/Users/FGreifeneder/Documents/tmp_processing_stuff/maps/tif/"

#search and filtering of radar data this will determine the sequence of files
dim_list <- list()
dimcntr <- 1
for (i in c(2005:2012)){

  paths <- list.files(paste(asar_path, "processing_", as.character(i), "/02_CONV/", sep=""), pattern="*.dim", full.names=TRUE)
  min_date <- as.numeric(as.Date(paste("05/01/", as.character(i), sep=""), format= "%m/%d/%Y")) + 2440588 #first of may
  max_date <- as.numeric(as.Date(paste("10/31/", as.character(i), sep=""), format= "%m/%d/%Y")) + 2440588 #thirtyfirst of october
  for (j in c(1:length(paths))){
    tmp <- read_NEST_DIM(paths[j])
    if (tmp$ACQUISITION_DATE >= min_date && tmp$ACQUISITION_DATE <= max_date){
      dim_list[[dimcntr]] <- tmp
      dimcntr <- dimcntr + 1
    } 
  }
}

dim_dates <- as.numeric()
for (i in c(1:length(dim_list))) dim_dates <- c(dim_dates, dim_list[[i]]$ACQUISITION_DATE)
dim_sort <- sort(dim_dates, index.return=TRUE)

GEOtopstack <- stack()
ASARstack <- stack()

#GEOtop simulation start- and end-data
minDate <- as.numeric(as.Date("05/01/2010", format= "%m/%d/%Y")) + 2440588
maxDate <- as.numeric(as.Date("06/01/2014", format= "%m/%d/%Y")) + 2440588

#Stacking ASAR and GEOtop data

for (i in c(1:length(dim_dates))){
  
  current_date <- dim_dates[i]
  current_gtopindex <- current_date-minDate+1
  geotoppath <- paste(model_path,"thetaliqL0001N", sprintf("%04d", current_gtopindex) ,".tif", sep="")
  fex <- file.exists(geotoppath)
  asarpath <- paste("C:/Users/FGreifeneder/Documents/tmp_processing_stuff/ASARPROCESSING/SVR_noNDVI_20141215/maps_extrapol/", dim_list[[dim_sort$ix[i]]]$DATASET_NAME, "_SMC.tif", sep="")
  asarfex <- file.exists(asarpath)
  
  #skip to next date if file doesn't exist
  if (fex == F | asarfex == F) next
  gtopfile <- raster(geotoppath)
  asarfile <- raster(asarpath)
  NAvalue(asarfile) <- -32768
  NAvalue(gtopfile) <- -9999
  
  asarlyr <- resample(asarfile, gtopfile, method="ngb", NAflag=-32768)
  gtoplyr <- resample(gtopfile, gtopfile, method="ngb", NAflag=-32768)
  
  GEOtopstack <- stack(GEOtopstack, gtoplyr)
  ASARstack <- stack(ASARstack, asarlyr)
                                         
  
}
                                             
                                             
#Analysing the data - regression parameters, correlation, etc.

corr <- matrix(-32768, nrow(ASARstack), ncol(ASARstack))
n <- matrix(-32768, nrow(ASARstack), ncol(ASARstack))
SVRerror <- matrix(-32768, nrow(ASARstack), ncol(ASARstack))
SVRparamlist <- list(x=as.numeric(), y=as.numeric(), SVR=NULL)

GEOtopstack_out <- GEOtopstack
values(GEOtopstack_out) <- -Inf

GEOtopvalid <- which(is.na(as.matrix(subset(GEOtopstack, 1))) == F, arr.ind=TRUE, useNames=FALSE)

nvalid <- length(GEOtopvalid)/2

for (ind in c(1:nvalid)){
  print(paste(ind, ":", nvalid, sep=""))
  GEOtopts <- extract(GEOtopstack, cellFromRowCol(GEOtopstack, GEOtopvalid[ind,1], GEOtopvalid[ind,2]))
  GEOtopts <- as.vector(GEOtopts)
  ASARts <- extract(ASARstack, cellFromRowCol(ASARstack, GEOtopvalid[ind,1], GEOtopvalid[ind,2]))
  ASARts <- as.vector(ASARts)
  val <- which(is.na(GEOtopts) == F & is.na(ASARts) == F, arr.ind=TRUE)
  if (length(val) > 1) {
#     SVRparams <- ksvm(GEOtopts, y=ASARts, 
#                       type="eps-svr", 
#                       kernel="rbfdot", 
#                       kpar="automatic",
#                       prob.model=TRUE,
#                       cross=10, 
#                       C=tuneResult$best.parameters["cost"],
#                       epsilon=tuneResult$best.parameters["epsilon"])
    SVRparams <- NULL
    tuneResult <- NULL
      
    SVRparams <- tryCatch({tmp <- cbind(ASARts[val], GEOtopts[val])
                           colnames(tmp) <- c('targ', 'feat')
                           tmp <- as.data.frame(tmp)
                           tuneResult <- tune(svm, targ ~ feat,  data = tmp, ranges = list(epsilon = seq(0,1,0.1), cost = 2^(2:9)))
                           rm(tmp)
                           print("Solving regression ...")
                           SVRparams <- ksvm(GEOtopts, y=ASARts, 
                                     type="eps-svr", 
                                     kernel="rbfdot", 
                                     kpar="automatic",
                                     prob.model=TRUE,
                                     cross=10,
                                     C=tuneResult$best.parameters["cost"],
                                     epsilon=tuneResult$best.parameters["epsilon"])}, 
                           warning = function(war){
                                        print(paste("One warning occured: ", war))
                           },
                           error = function(err){
                                        print(paste("An error occured:", err))
                                        f <- NULL
                                        return(f)
                                   }, 
                           finally={if (length(SVRparams)==0) next}
                          )
    #print(SVRparams)
    
    n[GEOtopvalid[ind,1], GEOtopvalid[ind,2]] <- length(val)
    SVRerror[GEOtopvalid[ind,1], GEOtopvalid[ind,2]] <- error(SVRparams)
    corr[GEOtopvalid[ind,1], GEOtopvalid[ind,2]] <- cor(GEOtopts, y=ASARts)
    
#     GEOtopts_out <- predict(SVRparams, GEOtopts)
#     GEOtopstack_out[cellFromRowCol(GEOtopstack, GEOtopvalid[ind,1], GEOtopvalid[ind,2])] <- GEOtopts_out
    SVRparamlist$x <- c(SVRparamlist$x, GEOtopvalid[ind,1])
    SVRparamlist$y <- c(SVRparamlist$y, GEOtopvalid[ind,2])
    SVRparamlist$SVR <- c(SVRparamlist$SVR, SVRparams)
#   } else {
#     GEOtopts_out <- GEOtopts
#     GEOtopts_out[] <- -Inf
#     GEOtopstack_out[cellFromRowCol(GEOtopstack, GEOtopvalid[ind,1], GEOtopvalid[ind,2])] <- GEOtopts_out
  }
}

      
      


