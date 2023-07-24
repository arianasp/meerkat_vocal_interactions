#CALL RATE VS SPEED

#Typical rate of calling as a function of...
# individual speed
# individual speed relative to group

library(lubridate)
library(zoo)
library(fields)
library(viridis)

#PARAMS
groupyears <- c('HM2017','HM2019','L2019')
dt <- 60 #time window
indir <- '/Users/Ari/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/'
outdir <- '/Users/Ari/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/plots_2023-07-20/callrate_dt60_R10_manualbins/'
min_n_tracked <- 5
n_dist_bins <- 8
max_dist <- 100 #distance above which we don't count dyadic distances for the mean dyadic distance change plots
R_neighbor < - 10 #local neighborhood radius
speed_bins <- c(0,1,2,5,10,25,50,400) #speed bins
neighbor_change_bins <- c(-10.5,-3.5,-2.5,-1.5,-.5,.5,1.5,2.5,3.5,10.5)

#LOOP OVER GROUPS AND STORE METRICS DATA
data <- list()

setwd(indir)
for(g in 1:length(groupyears)){
  
  groupyear <- groupyears[g]
  data[[g]] <- list()
  
  print(groupyear)
  
  #--------LOAD DATA----------
  load(paste0(groupyear,'_COORDINATES_all_sessions.RData'))
  calls <- read.csv(file = paste0('all_calls_sync_resolved_with_oor_2022-12-04.csv'), header=T, stringsAsFactors =F, sep = ',')
  load(paste0(groupyear,'_DATAPRESENCE_all_sessions.RData'))
  
  #get data associated with that groupyear - store with general names
  timeLine <- eval(as.name(paste(groupyear,'timeLine', sep = '_')))
  indInfo <- eval(as.name(paste(groupyear,'indInfo', sep = '_')))
  audioOn <- eval(as.name(paste(groupyear,'audioOn', sep = '_')))
  gpsOn <- eval(as.name(paste(groupyear,'gpsOn', sep = '_')))
  allX <- eval(as.name(paste(groupyear,'allX', sep = '_')))
  allY <- eval(as.name(paste(groupyear,'allY', sep = '_')))
  dayIdx <- eval(as.name(paste(groupyear,'dayIdx',sep='_')))
  
  #remove data when gps not on
  allX[!gpsOn] <- NA
  allY[!gpsOn] <- NA
  
  #get number of individuals and number of times
  n_inds <- nrow(allX)
  n_times <- ncol(allX)
  n_days <- length(dayIdx)-1
  
  #get number tracked
  n_tracked <- colSums(!is.na(allX))
  
  #-------------COMPUTE METRICS-----------

  #calculate individual speeds
  speeds <- matrix(NA, nrow = nrow(allX), ncol=ncol(allX))
  for(i in 1:n_inds){
    for(d in 1:n_days){
      t_idxs <- dayIdx[d]:(dayIdx[d+1]-1)
      t_fut <- t_idxs[(dt+1):length(t_idxs)]
      t_now <- t_idxs[1:(length(t_idxs)-dt)]
      dx <- allX[i,t_fut] - allX[i, t_now]
      dy <- allY[i,t_fut] - allY[i, t_now]
      speeds[i,t_now] <- sqrt(dx^2 + dy^2) / dt * 60
    }
  }
  
  #get dyadic distances
  dyad_dists <- array(NA, dim = c(n_inds,n_inds, n_times))
  for(i in 1:(n_inds-1)){
    for(j in i:n_inds){
      dyad_dists[i,j,] <- dyad_dists[j,i,] <- sqrt((allX[i,] - allX[j,])^2 + (allY[i,] - allY[j,])^2)
    }
    dyad_dists[i,i,] <- NA
  }
  dyad_dists[which(dyad_dists>max_dist)] <- NA #remove dyadic distances above the max dist threshold - don't include them in the analysis
  
  #get how many individuals in local neighborhood (within radius R_neighbor)
  n_neighbors <- apply(dyad_dists <= R_neighbor, c(1,3), sum, na.rm=T)
  n_neighbors[,which(n_tracked < min_n_tracked)] <- NA
  
  #calculate mean absolute change in dyadic distance relative to all neighbors within a max radius (mean_dyad_dist_change)
  #also calculate change in neighbors within R_neighbor (neighbor_change)
  mean_dyad_dist_change <- neighbor_change <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    for(d in 1:n_days){
      
      #get time indexes on that day
      t_idxs <- dayIdx[d]:(dayIdx[d+1]-1)
      
      #get future time indexes and current time indexes for computing changes
      t_fut <- t_idxs[(dt+1):length(t_idxs)]
      t_now <- t_idxs[1:(length(t_idxs)-dt)]
      
      #get mean change in dyadic distance relative to all neighbors < max_dist
      mean_dyad_dist_change[i, t_now] <- colMeans(abs(dyad_dists[i,,t_fut] - dyad_dists[i,,t_now]),na.rm=T)
      
      #get neighborhood change
      neighbor_change[i, t_now] <- n_neighbors[i, t_fut] - n_neighbors[i, t_now]
      
    }
  }
  
  #calculate number of calls in each second
  calls$datetime <- as.POSIXct(calls$t0GPS_UTC, tz = 'UTC')
  calls$tidx <- match(as.character(calls$datetime),as.character(as.POSIXct(timeLine,tz='UTC')))
  ccs <- calls[which(grepl('cc',calls$stn_call_type) & calls$pred_focalType=='F'),]
  sns <- calls[which(grepl('sn',calls$stn_call_type) & calls$pred_focalType=='F'),]
  
  cc_timeLine <- sn_timeLine <- matrix(0, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    
    cc_times <- ccs$tidx[which(ccs$ind == indInfo$code[i])]
    sn_times <- sns$tidx[which(sns$ind == indInfo$code[i])]
    cc_times <- cc_times[!is.na(cc_times)]
    sn_times <- sn_times[!is.na(sn_times)]
    for(t in 1:length(cc_times)){
      cc_timeLine[i,cc_times[t]] <- cc_timeLine[i,cc_times[t]]+1
    }
    for(t in 1:length(sn_times)){
      sn_timeLine[i,sn_times[t]] <- sn_timeLine[i,sn_times[t]]+1
    }
  }
  
  #remove data when audio was not labeled
  cc_timeLine[!audioOn] <- NA
  sn_timeLine[!audioOn] <- NA
  
  #get call rates in the time window dt
  cc_callrates <- sn_callrates <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    cc_callrates[i,] <- rollapply(cc_timeLine[i,], width = dt, FUN = sum, na.rm=T, fill=NA, align = 'left') / dt * 60
    sn_callrates[i,] <- rollapply(sn_timeLine[i,], width = dt, FUN = sum, na.rm=T, fill=NA, align = 'left') / dt * 60
  }
  
  #remove data for anything when audio was not on
  cc_callrates[!audioOn] <- NA
  sn_callrates[!audioOn] <- NA
  
  #normalize call probability by individual - use z score
  #TODO: How to normalize properly when data are very zero-inflated?
  cc_callrates_zscore <- sn_callrates_zscore <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    
    meancc <- mean(cc_callrates[i,], na.rm=T)
    sdcc <- sd(cc_callrates[i,],na.rm=T)
    cc_callrates_zscore[i,] <- (cc_callrates[i,] - meancc) / sdcc
    
    meansn <- mean(sn_callrates[i,], na.rm=T)
    sdsn <- sd(sn_callrates[i,],na.rm=T)
    sn_callrates_zscore[i,] <- (sn_callrates[i,] - meansn) / sdsn
    
  }
  
  data[[g]]$speeds <- speeds
  data[[g]]$mean_dyad_dist_change <- mean_dyad_dist_change
  data[[g]]$cc_callrates <- cc_callrates
  data[[g]]$sn_callrates <- sn_callrates
  data[[g]]$cc_callrates_zscore <- cc_callrates_zscore
  data[[g]]$sn_callrates_zscore <- sn_callrates_zscore
  data[[g]]$n_neighbors <- n_neighbors
  data[[g]]$neighbor_change <- neighbor_change
  data[[g]]$n_tracked <- n_tracked
  
}

#------------PLOTS------------
  
for(g in 1:length(groupyears)){
  
  groupyear <- groupyears[[g]]
  
  speeds <- data[[g]]$speeds 
  mean_dyad_dist_change <- data[[g]]$mean_dyad_dist_change 
  cc_callrates <- data[[g]]$cc_callrates 
  sn_callrates <- data[[g]]$sn_callrates 
  cc_callrates_zscore <- data[[g]]$cc_callrates_zscore 
  sn_callrates_zscore <- data[[g]]$sn_callrates_zscore 
  neighbor_change <- data[[g]]$neighbor_change
  n_neighbors <- data[[g]]$n_neighbors
  n_tracked <- data[[g]]$n_tracked
  
  #bins
  dyad_dist_change_bins <- quantile(mean_dyad_dist_change, seq(0,1,length.out=8),na.rm=T)
  n_neighbor_bins <- seq(-.5, max(n_neighbors,na.rm=T)+.5, 1)
  
  #----------Bivariate plots----
  #PLOT 1-2: Individual speed vs mean call rates
  quartz(height = 10, width = 6)
  par(mfrow=c(4,2))
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(speed_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    idxs <- which(speeds >= speed_bins[i] & speeds < speed_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (speed_bins[1:(length(speed_bins)-1)] + speed_bins[2:length(speed_bins)]) / 2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'CC call rate (calls / min)', log = 'x')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'SN call rate (calls / min)', log = 'x')
  
  #PLOT 3-4: Movement relative to group vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(dyad_dist_change_bins)-1)
  for(i in 1:(length(dyad_dist_change_bins)-1)){
    idxs <- which(mean_dyad_dist_change >= dyad_dist_change_bins[i] & mean_dyad_dist_change < dyad_dist_change_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (dyad_dist_change_bins[1:(length(dyad_dist_change_bins)-1)] + dyad_dist_change_bins[2:length(dyad_dist_change_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Mean dyad dist change (m)', ylab = 'Mean CC rate', log = 'x')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Mean dyad dist change (m)', ylab = 'Mean SN rate', log = 'x')
  
  #PLOT 5-6: Neighborhood change vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(neighbor_change_bins)-1)
  for(i in 1:(length(neighbor_change_bins)-1)){
    idxs <- which(neighbor_change >= neighbor_change_bins[i] & neighbor_change < neighbor_change_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (neighbor_change_bins[1:(length(neighbor_change_bins)-1)] + neighbor_change_bins[2:length(neighbor_change_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Neighbor change', ylab = 'Mean CC rate')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Neighbor change', ylab = 'Mean SN rate')

  #PLOT 7-8: N neighbors vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(n_neighbor_bins)-1)
  for(i in 1:(length(n_neighbor_bins)-1)){
    idxs <- which(n_neighbors >= n_neighbor_bins[i] & n_neighbors < n_neighbor_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (n_neighbor_bins[1:(length(n_neighbor_bins)-1)] + n_neighbor_bins[2:length(n_neighbor_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean CC rate')
  plot(mids, mean_sn_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean SN rate')
  
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_neighbors_',groupyear,'.pdf'))
  dev.off()
  
  #--------Bivariate plots - normalized -----
  #PLOT 1-2: Individual speed vs mean call rates - normed
  quartz(height = 10, width = 6)
  par(mfrow=c(4,2))
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(speed_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    idxs <- which(speeds >= speed_bins[i] & speeds < speed_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (speed_bins[1:(length(speed_bins)-1)] + speed_bins[2:length(speed_bins)]) / 2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean CC rate z-score', log = 'x')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean CC rate z-score', log = 'x')
  
  #PLOT 3-4: Movement relative to group vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(dyad_dist_change_bins)-1)
  for(i in 1:(length(dyad_dist_change_bins)-1)){
    idxs <- which(mean_dyad_dist_change >= dyad_dist_change_bins[i] & mean_dyad_dist_change < dyad_dist_change_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (dyad_dist_change_bins[1:(length(dyad_dist_change_bins)-1)] + dyad_dist_change_bins[2:length(dyad_dist_change_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Mean dyad dist change (m)', ylab = 'Mean CC rate z-score', log = 'x')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Mean dyad dist change (m)', ylab = 'Mean SN rate z-score', log = 'x')
  
  #PLOT 5-6: Neighborhood change vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(neighbor_change_bins)-1)
  for(i in 1:(length(neighbor_change_bins)-1)){
    idxs <- which(neighbor_change >= neighbor_change_bins[i] & neighbor_change < neighbor_change_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (neighbor_change_bins[1:(length(neighbor_change_bins)-1)] + neighbor_change_bins[2:length(neighbor_change_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = 'Neighbor change', ylab = 'Mean CC rate - zscore')
  plot(mids, mean_sn_rates, pch = 19, xlab = 'Neighbor change', ylab = 'Mean SN rate - zscore')
  
  #PLOT 7-8: N neighbors vs mean call rates
  mean_cc_rates <- mean_sn_rates <- rep(NA, length(n_neighbor_bins)-1)
  for(i in 1:(length(n_neighbor_bins)-1)){
    idxs <- which(n_neighbors >= n_neighbor_bins[i] & n_neighbors < n_neighbor_bins[i+1])
    mean_cc_rates[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
    mean_sn_rates[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
  }
  
  #make plot
  mids <- (n_neighbor_bins[1:(length(n_neighbor_bins)-1)] + n_neighbor_bins[2:length(n_neighbor_bins)])/2
  plot(mids, mean_cc_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean CC rate - zscore')
  plot(mids, mean_sn_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean SN rate - zscore')
  
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_neighbors_zscore',groupyear,'.pdf'))
  dev.off()
  
  #---------------3d plots---------------
  #3D plot - call rate vs speed and mean dyadic distance change
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(dyad_dist_change_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(dyad_dist_change_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      mean_dyad_dist_change >= dyad_dist_change_bins[j] &
                      mean_dyad_dist_change < dyad_dist_change_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
      }
    }
  }
  quartz()
  par(mfrow=c(2,2))
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Mean dyad dist change',main = 'CC')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Mean dyad dist change',main='SN')
  
  #3D plot - call rate vs speed and mean dyadic distance change - zscore
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(dyad_dist_change_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(dyad_dist_change_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      mean_dyad_dist_change >= dyad_dist_change_bins[j] &
                      mean_dyad_dist_change < dyad_dist_change_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
      }
    }
  }
  
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Mean change in dyad dist',main = 'CC - zscore')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Mean change in dyad dist',main='SN - zscore')

  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_rel_dyad_dist_change_2d_',groupyear,'.pdf'))
  dev.off()
  
  #---------------3d plots neighbor change---------------
  #3D plot - call rate vs speed and mean dyadic distance change
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(neighbor_change_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(neighbor_change_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      neighbor_change >= neighbor_change_bins[j] &
                      neighbor_change < neighbor_change_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
      }
    }
  }
  quartz()
  par(mfrow=c(2,2))
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Neighbor change',main = 'CC')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Neighbor change',main='SN')
  
  #3D plot - call rate vs speed and mean dyadic distance change - zscore
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(neighbor_change_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(neighbor_change_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      neighbor_change >= neighbor_change_bins[j] &
                      neighbor_change < neighbor_change_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
      }
    }
  }
  
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Neighbor change',main = 'CC - zscore')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='Neighbor change',main='SN - zscore')
  
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_neighbor_change_2d_',groupyear,'.pdf'))
  dev.off()
  
  #---------------3d plots number of neighbors---------------
  #3D plot - call rate vs speed and mean dyadic distance change
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(n_neighbor_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(n_neighbor_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      n_neighbors >= n_neighbor_bins[j] &
                      n_neighbors < n_neighbor_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
      }
    }
  }
  quartz()
  par(mfrow=c(2,2))
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main = 'CC')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main='SN')
  
  #3D plot - call rate vs speed and mean dyadic distance change - zscore
  cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(n_neighbor_bins)-1)
  for(i in 1:(length(speed_bins)-1)){
    for(j in 1:(length(n_neighbor_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] &
                      speeds < speed_bins[i+1] &
                      n_neighbors >= n_neighbor_bins[j] &
                      n_neighbors < n_neighbor_bins[j+1])
      nsamps[i,j] <- length(idxs)
      if(length(idxs)>1000){
        cc_rates_2d[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
        sn_rates_2d[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
      }
    }
  }
  
  image.plot(z=cc_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main = 'CC - zscore')
  image.plot(z=sn_rates_2d, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main='SN - zscore')
  
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_n_neighbors_2d_',groupyear,'.pdf'))
  dev.off()
  
  
}

#Look at speed results when at least half group giving close calls (filter out non-foraging periods)
#Call rate vs. number of neighbors and surroundedness (circular variance)
#Call rate vs number of neighbors and their call rates (?)
