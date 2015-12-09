#returns a data frame, containnig the scaled pixel data (based on the SMAP grid)

#smap pixels:
#
#1: B2
#2: (B3)
#3: B1
#4: M1
#5: M5
#6:
#7:
#8: M3/M4
#9: M2
#10:

library("insol")

extrScaledSMC <- function(path){
  
  scaled_SMC <- list()
  
  for (i in 1:10){
    
    scaled_SMC[[i]] <- read.table(paste(path,"/vwc_scaled_SMAP",sprintf("%02d", i),".txt", sep=""),
                                  sep=",", header=T, stringsAsFactors=F, quote="")
    
  }
  
  startDay <- min(min(scaled_SMC[[1]][,1]), 
                  min(scaled_SMC[[2]][,1]),
                  min(scaled_SMC[[3]][,1]),
                  min(scaled_SMC[[4]][,1]),
                  min(scaled_SMC[[5]][,1]),
                  min(scaled_SMC[[6]][,1]),
                  min(scaled_SMC[[7]][,1]),
                  min(scaled_SMC[[8]][,1]),
                  min(scaled_SMC[[9]][,1]),
                  min(scaled_SMC[[10]][,1]))
  
  endDay <- max(max(scaled_SMC[[1]][,1]), 
                max(scaled_SMC[[2]][,1]),
                max(scaled_SMC[[3]][,1]),
                max(scaled_SMC[[4]][,1]),
                max(scaled_SMC[[5]][,1]),
                max(scaled_SMC[[6]][,1]),
                max(scaled_SMC[[7]][,1]),
                max(scaled_SMC[[8]][,1]),
                max(scaled_SMC[[9]][,1]),
                max(scaled_SMC[[10]][,1]))
  
  df_scaled <- data.frame(date=startDay:endDay, 
                          p1=rep(-1,endDay-startDay+1),
                          p2=rep(-1,endDay-startDay+1),
                          p3=rep(-1,endDay-startDay+1),
                          p4=rep(-1,endDay-startDay+1),
                          p5=rep(-1,endDay-startDay+1),
                          p6=rep(-1,endDay-startDay+1),
                          p7=rep(-1,endDay-startDay+1),
                          p8=rep(-1,endDay-startDay+1),
                          p9=rep(-1,endDay-startDay+1),
                          p10=rep(-1,endDay-startDay+1))
  
  for (i in 1:(endDay-startDay)){
    
    if (any(scaled_SMC[[1]][,1] == df_scaled$date[i])) { 
      df_scaled$p1[i] <- scaled_SMC[[1]][which(scaled_SMC[[1]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[2]][,1] == df_scaled$date[i])) { 
      df_scaled$p2[i] <- scaled_SMC[[2]][which(scaled_SMC[[2]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[3]][,1] == df_scaled$date[i])) { 
      df_scaled$p3[i] <- scaled_SMC[[3]][which(scaled_SMC[[3]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[4]][,1] == df_scaled$date[i])) { 
      df_scaled$p4[i] <- scaled_SMC[[4]][which(scaled_SMC[[4]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[5]][,1] == df_scaled$date[i])) { 
      df_scaled$p5[i] <- scaled_SMC[[5]][which(scaled_SMC[[5]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[6]][,1] == df_scaled$date[i])) { 
      df_scaled$p6[i] <- scaled_SMC[[6]][which(scaled_SMC[[6]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[7]][,1] == df_scaled$date[i])) { 
      df_scaled$p7[i] <- scaled_SMC[[7]][which(scaled_SMC[[7]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[8]][,1] == df_scaled$date[i])) { 
      df_scaled$p8[i] <- scaled_SMC[[8]][which(scaled_SMC[[8]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[9]][,1] == df_scaled$date[i])) { 
      df_scaled$p9[i] <- scaled_SMC[[9]][which(scaled_SMC[[9]][,1] == df_scaled$date[i]),2] 
    }
    if (any(scaled_SMC[[10]][,1] == df_scaled$date[i])) { 
      df_scaled$p10[i] <- scaled_SMC[[10]][which(scaled_SMC[[10]][,1] == df_scaled$date[i]),2] 
    }
    
  }
  
  df_scaled$date <- as.Date(JD(df_scaled$date, inverse=T))
  
  return(df_scaled)
  
}