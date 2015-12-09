#create daily averages from SMC time series

dailySMCavg <- function(path){
  
  load(path)
  ts30 <- combinedTable
  val <- which(ts30$SM05_1 != -7777 | ts30$SM05_2 != -7777)
  ts30 <- ts30[val,]
  ts30$SM05_1[ts30$SM05_1 == -7777] <- NA
  ts30$SM05_2[ts30$SM05_2 == -7777] <- NA
  
  dailySMC <- data.frame(date=as.character(), mean_SM05=as.numeric(), stringsAsFactors = F)
  
  days <- as.Date(paste(ts30$Yr,ts30$Mo,ts30$Dy,sep="-"))
  uniqueDays <- unique(days)
  cntr <- 1
  
  for (cDay in uniqueDays){
    subID <- which(days == cDay)
    dailySMC[cntr,1] <- as.character.Date(uniqueDays[cntr])
    dailySMC[cntr,2] <- mean(rowMeans(rbind(ts30$SM05_1[subID],ts30$SM05_2[subID]), na.rm=T), na.rm=T)
    cntr <- cntr + 1
    
  }
  
  #names(dailySMC) <- c("Date", "mean_SM05")
  
  return(dailySMC)
  
  
}

dailySMCavg_scaled <- function(path){
  
  load(path)
  ts30 <- combinedTable
  val <- which(ts30$WASM_1 != -7777 | ts30$WASM_2 != -7777)
  ts30 <- ts30[val,]
  ts30$WASM_1[ts30$WASM_1 == -7777] <- NA
  ts30$WASM_2[ts30$WASM_2 == -7777] <- NA
  
  dailySMC <- data.frame(date=as.character(), mean_SM05=as.numeric(), stringsAsFactors = F)
  
  days <- as.Date(paste(ts30$Yr,ts30$Mo,ts30$Dy,sep="-"))
  uniqueDays <- unique(days)
  cntr <- 1
  
  for (cDay in uniqueDays){
    subID <- which(days == cDay)
    dailySMC[cntr,1] <- as.character.Date(uniqueDays[cntr])
    dailySMC[cntr,2] <- mean(rowMeans(rbind(ts30$WASM_1[subID],ts30$WASM_2[subID]), na.rm=T), na.rm=T)
    cntr <- cntr + 1
    
  }
  
  #names(dailySMC) <- c("Date", "mean_SM05")
  
  return(dailySMC)
  
  
}