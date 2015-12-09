#this function computes the mutual information index between one column in a matrix an all other

library(entropy)

computeMI <- function(data, targetIndex){

	target <- data[,targetIndex]
	featureIndex <- 1:ncol(data)
	#featureIndex <- featureIndex[-targetIndex]
	#print(featureIndex)

	#discretize target
	#target <- (target - min(target)) / (max(target) - min(target))
	target <- freqs.empirical(discretize(target, numBins=50))

	MI <- numeric()

	for (i in featureIndex){
		#discretize vectors
		#feature <- data[,i]
		#feature <- (feature - min(feature)) / (max(feature) - min(feature))
		feature <- freqs.empirical(discretize(data[,i], numBins=50))
		MI <- c(MI, mi.plugin(rbind(target, feature)))
	}
	
	featureNames <- colnames(data)
	#print(featureNames)
	#featureNames <- featureNames[-targetIndex]
	#print(featureNames)
	#print(MI)
	names(MI) <- featureNames
	return(MI)
}

data20m <- read.csv("X:/Workspaces/GrF/01_Data/RADARSAT2/quadpol2014/parameters/20140715_reprocessing/parameters_20140715_CSV_CLSTRD.csv")
data20m <- data20m[,4:82]

output <- as.data.frame(computeMI(data20m, 1))
for (i in c(2:length(data20m))){
	output <- cbind(output, computeMI(data20m, i))
}


write.csv(output, "X:/Workspaces/GrF/Processing/RADARSAT/cp_reprocessing/july_nest_20m/MI_SMC_features_CLSTRD_Rout.csv")
