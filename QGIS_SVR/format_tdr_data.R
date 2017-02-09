# this scripts useses the output of the TDR, as a csv file, an averages based on zones

#INPUT
table <- read.csv("Q:\\ESA_TIGER\\processed\\in-situ\\field_campaign_17042016\\tdrmerged_18042016.csv")

table$Moisture[which(table$Moisture > 100)] <- NA
table_avg <- aggregate(cbind(Moisture,Latitude,Longitude)~Zone, data = table, FUN=mean)
#table_avg <- aggregate(cbind(Moisture, La)~ID, data = table, FUN=mean)

#OUTPUT
write.csv(table_avg, file="Q:\\ESA_TIGER\\processed\\in-situ\\field_campaign_17042016\\tdrmerged_18042016_aggregated.csv", quote=F)
