#Time differences between pairs of calls of individuals within range dist.thresh

#directory where gps data is stored for the project
#gps.datadir <- '/Volumes/EAS_shared/meerkat/working/processed/movement/' #SERVER
gps.datadir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/' #LOCAL

#direcotry where audio labeling data is stored for the project
#audio.datadir <- '/Volumes/EAS_shared/meerkat/working/processed/acoustic' #SERVER
audio.datadir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/' #LOCAL

#directory where code is stored for this project
codedir <- '~/Dropbox/code_ari/meerkat_vocal_interactions'

#directory of where to save results
savedir <- '~/Dropbox/meerkats/results/call_interactions/'

#filename of general meerkat functions
general.funcs.filename <- 'meerkat_functions.R'

#maximum time lag to consider (time since another individual called)
max.lag <- 10

#maximum distance to consider
max.dist <- 10

calltype <- 'cc'

#list of group years
groupyears <- c('HM2017', 'HM2019', 'L2019')

#file names
#audio.file <- 'full_labelfile_conflicts_resolved.csv'
gps.files <- paste(groupyears, 'COORDINATES_all_sessions.RData', sep = '_')

#audio.file <- '/Volumes/EAS_shared/meerkat/working/processed/acoustic/resolve_conflicts/all_calls_sync_resolved_with_oor_2022-12-04.csv' #SERVER
audio.file <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/all_calls_sync_resolved_with_oor_2022-12-04.csv' #LOCAL

#------------------LOAD DATA-----------------------

timestamp()
print('loading data')

#load audio data and datapresence data
setwd(audio.datadir)
calls.all <- read.csv(audio.file, header=T, sep=',', stringsAsFactors=F)
for(g in 1:length(groupyears)){
  load(paste(groupyears[g],'DATAPRESENCE_all_sessions.RData',sep='_'))
}

setwd(gps.datadir)
#load gps data
for(i in 1:length(gps.files)){
  load(gps.files[i])
}

#-------------------PROCESS-----------------------
#convert times to POSIXlt and then numeric
calls.all$t0 <- as.numeric(as.POSIXlt(calls.all$t0GPS_UTC, tz = 'UTC'))
calls.all$tf <- as.numeric(as.POSIXlt(calls.all$tendGPS_UTC, tz = 'UTC'))

#get location of calls


#get time differences between nearby individuals (using times of emission of call for the location of each ind)
time.diffs <- c()
dists <- c()

dates <- unique(calls.all$date)

for(d in 1:length(dates)){
  
  timestamp()
  print(d)
  date <- dates[d]
  calls.date <- calls.all[which(calls.all$date==date),]
  
  if(nrow(calls.date)==0){
    next
  }
  
  for(i in 1:nrow(calls.date)){
    xi <- calls.date$x_emitted[i]
    yi <- calls.date$y_emitted[i]
    t0 <- calls.date$t0[i]
    
    dist.to.caller <- sqrt((calls.date$x_emitted - xi)^2 + (calls.date$y_emitted - yi)^2)
    time.to.call <- calls.date$t0 - t0
    
    idxs <- which(dist.to.caller <= max.dist & time.to.call <= max.lag & time.to.call > 0 & dist.to.caller > 0)
    
    if(length(idxs)==0){
      next
    }
    
    time.diffs <- c(time.diffs, time.to.call[idxs])
    dists <- c(dists, dist.to.caller[idxs])
    
  }
  
  
}




