csvs <- unique(calls.all$csvFileName)

missing_stop <- c()
missing_start <- c()

for(i in 1:length(csvs)){
  tmp <- calls.all[which(calls.all$csvFileName == csvs[i]),]
  starts <- sum(tmp$callType == 'start')
  if(starts == 0){
    missing_start <- c(missing_start, csvs[i])
  }
  stops <- sum(tmp$callType == 'stop')
  if(stops == 0){
    missing_stop <- c(missing_stop, csvs[i])
  }
  
}

missing_start <- as.data.frame(missing_start)
write.table(missing_start, file = '/Users/Ari/Desktop/missing_start.csv',sep=',')
missing_end <- as.data.frame(missing_stop)
write.table(missing_stop, file = '/Users/Ari/Desktop/missing_stop.csv',sep=',')
