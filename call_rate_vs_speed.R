#CALL RATE VS SPEED

#Typical rate of calling as a function of...
# individual speed
# individual speed relative to group

library(lubridate)
library(zoo)
library(fields)
library(viridis)

#FUNCTIONS
source('~/Dropbox/code_ari/Meerkat_leadership_hackathon/INFLUENCE_PAPER/scripts/functions.R')

#PARAMS
groupyears <- c('HM2017','HM2019','L2019')
dt <- 10 #time window
indir <- '/Users/Ari/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/'
outdir <- '/Users/Ari/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/plots_2023-07-20/callrate_dt10_manualbins/'
min_n_tracked <- 5
n_speed_bins <- 9
n_dist_bins <- 8
max_dist <- 100 #which distance to use above which we don't count dyadic distances


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
  
  #-------------COMPUTE METRICS-----------
  
  #get group centroid
  centrX <- colMeans(allX, na.rm=T)
  centrY <- colMeans(allY, na.rm=T)
  
  #remove data for fewer than min_n_tracked individuals tracked
  n_tracked <- colSums(!is.na(allX))
  centrX[which(n_tracked < min_n_tracked)] <- NA
  centrY[which(n_tracked < min_n_tracked)] <- NA
  
  #get spatially discretized headings and relative positions etc
  relPos <- relativePos(x = allX, y = allY, discrSpatial = T, step = 20, futur=F)
  
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
  }
  dyad_dists[which(dyad_dists>max_dist)] <- NA #remove dyadic distances above the max dist threshold - don't include them in the analysis
  
  #calculate movement relative to group centroid
  #change in distance from centroid (dist_change) 
  #change in relative position relative to moving reference frame of group (rel_pos_change)
  rel_pos_change <- dist_change <- dist_from_centr <- mean_dyad_dist_change <- matrix(NA, nrow = n_inds, ncol = n_times)
  for(i in 1:n_inds){
    for(d in 1:n_days){
      
      #get time indexes on that day
      t_idxs <- dayIdx[d]:(dayIdx[d+1]-1)
      
      #get future time indexes and current time indexes for computing changes
      t_fut <- t_idxs[(dt+1):length(t_idxs)]
      t_now <- t_idxs[1:(length(t_idxs)-dt)]
      
      #distance of the individual from the centroid over time
      dist_from_centr[i,t_idxs] <- sqrt((allX[i,t_idxs] - centrX[t_idxs])^2 + (allY[i,t_idxs] - centrY[t_idxs])^2)
      
      #change in distance from centroid
      dist_change[i, t_now] <- (dist_from_centr[i,t_fut] - dist_from_centr[i,t_now]) / dt * 60
      
      #change in relative position
      rel_pos_change[i, t_now] <- sqrt((relPos$relative_ind_X[i, t_fut] - relPos$relative_ind_X[i, t_now])^2 + (relPos$relative_ind_Y[i,t_fut] - relPos$relative_ind_Y[i,t_now])^2)
      
      #get mean change in dyadic distance relative to all neighbors < max_dist
      mean_dyad_dist_change[i, t_now] <- colMeans(abs(dyad_dists[i,,t_fut] - dyad_dists[i,,t_now]),na.rm=T)
      
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
  
}
  
#----------------GET BINS--------------------
#log_speed_bins <- seq(quantile(log_speeds,0.01,na.rm=T), quantile(log_speeds,0.99,na.rm=T), length.out = n_speed_bins)
#dist_centr_bins <- seq(0, quantile(dist_from_centr, 0.99, na.rm=T), length.out = n_dist_bins)
#dist_change_bins <- seq(quantile(dist_change, 0.01,na.rm=T), quantile(dist_change, 0.99, na.rm=T), length.out=n_dist_bins)
#relX_bins <- c(-100,-50,-25,-10, -5, 0, 5, 10, 25, 50, 100)
#rel_pos_change_bins <- seq(quantile(log(rel_pos_change), 0.01, na.rm=T), quantile(log(rel_pos_change), 0.99, na.rm=T), length.out = n_dist_bins)
#speed_bins <- quantile(speeds, seq(0,1,length.out=n_speed_bins),na.rm=T)
speed_bins <- c(0,1,2,5,10,25,50,400)

#------------PLOTS------------
  
for(g in 1:length(groupyears)){
  
  groupyear <- groupyears[[g]]
  
  speeds <- data[[g]]$speeds 
  mean_dyad_dist_change <- data[[g]]$mean_dyad_dist_change 
  cc_callrates <- data[[g]]$cc_callrates 
  sn_callrates <- data[[g]]$sn_callrates 
  cc_callrates_zscore <- data[[g]]$cc_callrates_zscore 
  sn_callrates_zscore <- data[[g]]$sn_callrates_zscore 
  
  dyad_dist_change_bins <- quantile(mean_dyad_dist_change, seq(0,1,length.out=8),na.rm=T)
  
  #----------Bivariate plots----
  #PLOT 1-2: Individual speed vs mean call rates
  quartz()
  par(mfrow=c(2,2))
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
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_rel_dyad_dist_change_',groupyear,'.pdf'))
  dev.off()
  
  #--------Bivariate plots - normalized -----
  #PLOT 1-2: Individual speed vs mean call rates - normed
  quartz()
  par(mfrow=c(2,2))
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
  dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_and_rel_dyad_dist_change_zscore',groupyear,'.pdf'))
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
}

# #PLOT: Individual speed vs mean call rates - by individual
# mean_cc_rates <- mean_sn_rates <- nsamp <- matrix(NA, nrow = n_inds, ncol = length(log_speed_bins)-1)
# for(i in 1:(length(log_speed_bins)-1)){
#   for(ind in 1:n_inds){
#     idxs <- which(log_speeds[ind,] >= log_speed_bins[i] & log_speeds[ind,] < log_speed_bins[i+1])
#     mean_cc_rates[ind, i] <- mean(cc_callrates[ind,idxs]>0, na.rm=T)
#     mean_sn_rates[ind, i] <- mean(sn_callrates[ind,idxs]>0, na.rm=T)
#     nsamp[ind,i] <- length(idxs)
#   }
# }
# quartz()
# par(mfrow=c(1,2))
# plot(NULL, xlim = c(min(log_speed_bins),max(log_speed_bins)), ylim = c(0, max(mean_cc_rates,na.rm=T)), xlab = 'Log speed (m/min)', ylab = 'P(CC)')
# for(ind in 1:n_inds){
#   lines(log_speed_bins[2:length(log_speed_bins)], mean_cc_rates[ind,], col = 'gray')
# }
# 
# plot(NULL, xlim = c(min(log_speed_bins),max(log_speed_bins)), ylim = c(0, max(mean_sn_rates,na.rm=T)), xlab = 'Log speed (m/min)', ylab = 'P(SN)')
# for(ind in 1:n_inds){
#   lines(log_speed_bins[2:length(log_speed_bins)], mean_sn_rates[ind,], col = 'gray')
# }
# 
# # #PLOT: Front-back position vs call rate - aggregated across individuals and by individual
# # mean_cc_rates_all <- mean_sn_rates_all <- nsamp_all <- rep(NA, length(relX_bins)-1)
# # for(i in 1:length(mean_cc_rates_all)){
# #   idxs <- which(relPos$relative_ind_X >= relX_bins[i] & relPos$relative_ind_X < relX_bins[i+1])
# #   nsamp[i] <- length(idxs)
# #   mean_cc_rates_all[i] <- mean(cc_callrates[idxs]>0,na.rm=T)
# #   mean_sn_rates_all[i] <- mean(sn_callrates[idxs]>0,na.rm=T)
# # }
# # 
# # mean_cc_rates <- mean_sn_rates <- nsamp <- matrix(NA, nrow = n_inds, ncol = length(relX_bins)-1)
# # for(i in 1:n_inds){
# #   for(j in 1:length(mean_cc_rates_all)){
# #     idxs <- which(relPos$relative_ind_X[i,] >= relX_bins[j] & relPos$relative_ind_X[i,] < relX_bins[j+1])
# #     nsamp[i,j] <- length(idxs)
# #     if(length(idxs)>0){
# #       mean_cc_rates[i,j] <- mean(cc_callrates[i,idxs]>0,na.rm=T)
# #       mean_sn_rates[i,j] <- mean(sn_callrates[i,idxs]>0,na.rm=T)
# #     }
# #   }
# # }
# # 
# # quartz()
# # par(mfrow=c(2,2))
# # mids <- (relX_bins[1:(length(relX_bins)-1)] + relX_bins[2:length(relX_bins)])/2
# # plot(mids, mean_cc_rates_all, pch = 19, xlab = 'Front-back position (m)', ylab = 'P(CC)')
# # plot(mids, mean_sn_rates_all, pch = 19, xlab = 'Front-back position (m)', ylab = 'P(SN)')
# # 
# # plot(NULL, xlim = c(min(relX_bins), max(relX_bins)),ylim = c(min(mean_cc_rates,na.rm=T), max(mean_cc_rates,na.rm=T)),xlab = 'Front-back position (m)', ylab = 'P(CC)')
# # for(i in 1:n_inds){
# #   lines(mids, mean_cc_rates[i,],lwd=2,col='gray')
# # }
# # 
# # plot(NULL, xlim = c(min(relX_bins), max(relX_bins)),ylim = c(min(mean_sn_rates,na.rm=T), max(mean_sn_rates,na.rm=T)),xlab = 'Front-back position (m)', ylab = 'P(SN)')
# # for(i in 1:n_inds){
# #   lines(mids, mean_sn_rates[i,],lwd=2,col='gray')
# # }
# # 
# 
# #---------------------multivariable heat map plots-------------
# #PLOT 4: 3D plot - call rate vs dist from centroid change and speed
# cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(log_speed_bins)-1, ncol = length(dist_change_bins)-1)
# for(i in 1:(length(log_speed_bins)-1)){
#   for(j in 1:(length(dist_change_bins)-1)){
#     idxs <- which(log_speeds >= log_speed_bins[i] &
#                     log_speeds < log_speed_bins[i+1] &
#                     dist_change >= dist_change_bins[j] &
#                     dist_change < dist_change_bins[j+1])
#     nsamps[i,j] <- length(idxs)
#     if(length(idxs)>1000){
#       cc_rates_2d[i,j] <- mean(cc_callrates[idxs]>0,na.rm=T)
#       sn_rates_2d[i,j] <- mean(sn_callrates[idxs]>0,na.rm=T)
#     }
#   }
# }
# quartz()
# par(mfrow=c(1,2))
# image.plot(log_speed_bins, dist_change_bins, cc_rates_2d, zlim = c(0, max(cc_rates_2d,na.rm=T)), col = viridis(256), xlab = 'Log speed', ylab ='Change in dist from centroid',main = 'CC')
# image.plot(log_speed_bins, dist_change_bins, sn_rates_2d, zlim = c(0, max(cc_rates_2d,na.rm=T)), col = viridis(256), xlab = 'Log speed', ylab ='Change in dist from centroid',main='SN')
# 
# #PLOT 5: 3D plot - call rate vs dist from centroid and speed
# dist_centr_bins <- seq(0, quantile(dist_from_centr, 0.95, na.rm=T), length.out = 12)
# cc_rates_2d <- sn_rates_2d <- nsamps <- matrix(NA, nrow = length(log_speed_bins)-1, ncol = length(dist_centr_bins)-1)
# for(i in 1:(length(log_speed_bins)-1)){
#   for(j in 1:(length(dist_centr_bins)-1)){
#     idxs <- which(log_speeds >= log_speed_bins[i] &
#                     log_speeds < log_speed_bins[i+1] &
#                     dist_from_centr >= dist_centr_bins[j] &
#                     dist_from_centr < dist_centr_bins[j+1])
#     nsamps[i,j] <- length(idxs)
#     if(length(idxs)>1000){
#       cc_rates_2d[i,j] <- mean(cc_callrates[idxs]>0,na.rm=T)
#       sn_rates_2d[i,j] <- mean(sn_callrates[idxs]>0,na.rm=T)
#     }
#   }
# }
# quartz()
# par(mfrow=c(1,2))
# image.plot(log_speed_bins, dist_centr_bins, cc_rates_2d, zlim = c(0, max(cc_rates_2d,na.rm=T)), col = viridis(256), xlab = 'Log speed', ylab ='Dist from centroid',main = 'CC')
# image.plot(log_speed_bins, dist_centr_bins, sn_rates_2d, zlim = c(0, max(cc_rates_2d,na.rm=T)), col = viridis(256), xlab = 'Log speed', ylab ='Dist from centroid',main='SN')
# 

