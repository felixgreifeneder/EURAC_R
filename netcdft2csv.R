library(ncdf4)

#file_list <- list.files("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/", pattern="*.nc", recursive=T)
file_list <- "P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/roi_5/modis_vi.nc"

for (findex in c(1:length(file_list))){
	output_list <- list()
	#datafile <- try(nc_open(paste("P:/01_Projects/EUMETSAT_SCA Cross-pol/05_datasets/datasets_v2/", file_list[findex], sep="")))
  datafile <- try(nc_open(file_list))
	if(class(datafile) == "try-error") {
    print("Error reading: ")
    print(file_list[findex])
    next
	}
	datanames <- names(datafile$var)
	for (nindex in c(1:length(datanames))){
		tmp <- ncvar_get(datafile, varid=datanames[nindex])
		output_list <- c(output_list, list(tmp))
	}

	names(output_list) <- datanames
	maxlen <-0

	for (lind in output_list){
		if (maxlen < length(lind)) {
			maxlen <- length(lind)
		}
	}

	output_matrix <- mat.or.vec(maxlen, length(datanames))

	for (i in c(1:length(output_list))){
		if (length(output_list[[i]]) == 1) {
			output_matrix[,i] <- rep(output_list[[i]], maxlen)
		} else {
			output_matrix[,i] <- output_list[[i]]
		}
	}
	
  #	output_matrix <- as.character(output_matrix)
	) 
							(unlist(strsplit(file_list[findex], split=".", fixed=T))[1]), ".csv", sep=""), row.names=F, col.names=datanames, sep=",", quote=FALSE)
  nc_close(datafile)
}