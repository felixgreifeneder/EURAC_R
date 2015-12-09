#do parallellised prediction
library(foreach)
library(doSNOW)


parallel_predictions<-function(fit,testing)
{
  cl<-makeCluster(4)
  registerDoSNOW(cl)
  num_splits<-4
  split_testing<-sort(rank(1:nrow(testing))%%4)
  predictions<-foreach(i=unique(split_testing),
                       .combine=c,.packages=c("e1071")) %dopar% {
                         as.numeric(predict(fit,newdata=testing[split_testing==i,]))
                       }
  stopCluster(cl)
  return(predictions)
}