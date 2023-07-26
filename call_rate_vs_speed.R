#CALL RATE VS SPEED

#Typical rate of calling as a function of...
# individual speed
# number of neighbors
# surroundedness 

library(lubridate)
library(zoo)
library(fields)
library(viridis)

#FLAGS
compute <- T
make_plots <- F


if(compute){
  #PARAMS
  groupyears <- c('HM2017','HM2019','L2019')
  dt <- 10 #time window
  indir <- '/Users/Ari/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/'
  outdir <- '/Users/Ari/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/plots_2023-07-26/callrate_dt10_dtforage60_R10/'
  min_n_tracked <- 5
  R_neighbor <- 10 #local neighborhood radius - for computing number of neighbors
  dt_forage <- 60 #time window to use to specfiy whether the group is foraging (if at least half of the group gives a cc during this window)
  
  #bins
  speed_bins <- c(0,1,2,5,10,25,50,100,400) #speed bins
  surroundedness_bins <- seq(0,1,length.out = 10) 
  neighbor_change_bins <- c(-10.5, -3.5, -2.5, -1.5, 0.5, 1.5, 2.5, 3.5, 10.5)
  local_call_bins <- c(0, 1, 2, 5, 10, 20, 50, 200)
  
  #------LOOP OVER GROUPS AND STORE METRICS DATA------
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
  
    #get how many individuals in local neighborhood (within radius R_neighbor)
    n_neighbors <- apply(dyad_dists <= R_neighbor, c(1,3), sum, na.rm=T)
    n_neighbors[,which(n_tracked < min_n_tracked)] <- NA
    
    #get surroundedness - circular variance of vectors pointing to all neighbors
    surroundedness <- matrix(NA, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      dx <- allX - matrix(rep(allX[i,], each = n_inds), nrow = n_inds, ncol = n_times)
      dy <- allY - matrix(rep(allY[i,], each = n_inds), nrow = n_inds, ncol = n_times)
      dx[i,] <- NA
      dy[i,] <- NA
      norms <- sqrt(dx^2 + dy^2)
      dx_u <- dx / norms
      dy_u <- dy / norms
      surroundedness[i,] <- -sqrt(colSums(dx_u, na.rm=T)^2 + colSums(dy_u, na.rm=T)^2) / colSums(!is.na(dx_u)) + 1
    }
    
    #calculate change in neighbors within R_neighbor (neighbor_change)
    neighbor_change <- matrix(NA, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      for(d in 1:n_days){
        
        #get time indexes on that day
        t_idxs <- dayIdx[d]:(dayIdx[d+1]-1)
        
        #get future time indexes and current time indexes for computing changes
        t_fut <- t_idxs[(dt+1):length(t_idxs)]
        t_now <- t_idxs[1:(length(t_idxs)-dt)]
  
        #get neighborhood change
        neighbor_change[i, t_now] <- n_neighbors[i, t_fut] - n_neighbors[i, t_now]
        
      }
    }
    
    #calculate number of calls in each second
    calls$datetime <- as.POSIXct(calls$t0GPS_UTC, tz = 'UTC')
    calls$tidx <- match(as.character(calls$datetime),as.character(as.POSIXct(timeLine,tz='UTC')))
    ccs <- calls[which(grepl('cc',calls$stn_call_type) & calls$pred_focalType=='F'),]
    sns <- calls[which(grepl('sn',calls$stn_call_type) & calls$pred_focalType=='F'),]
    als <- calls[which(grepl('al',calls$stn_call_type) & calls$pred_focalType=='F'),]
    
    cc_timeLine <- sn_timeLine <- al_timeLine <- matrix(0, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      
      cc_times <- ccs$tidx[which(ccs$ind == indInfo$code[i])]
      sn_times <- sns$tidx[which(sns$ind == indInfo$code[i])]
      al_times <- als$tidx[which(als$ind == indInfo$code[i])]
      cc_times <- cc_times[!is.na(cc_times)]
      sn_times <- sn_times[!is.na(sn_times)]
      al_times <- al_times[!is.na(al_times)]
      for(t in 1:length(cc_times)){
        cc_timeLine[i,cc_times[t]] <- cc_timeLine[i,cc_times[t]]+1
      }
      for(t in 1:length(sn_times)){
        sn_timeLine[i,sn_times[t]] <- sn_timeLine[i,sn_times[t]]+1
      }
      for(t in 1:length(al_times)){
        al_timeLine[i,al_times[t]] <- al_timeLine[i,al_times[t]]+1
      }
    }
    
    #remove data when audio was not labeled
    cc_timeLine[!audioOn] <- NA
    sn_timeLine[!audioOn] <- NA
    al_timeLine[!audioOn] <- NA
    
    #get call rates in the time window dt
    cc_callrates <- sn_callrates <- cc_prev_callrates <- al_prev_callrates <- sn_prev_callrates <- matrix(NA, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      cc_callrates[i,] <- rollapply(cc_timeLine[i,], width = dt, FUN = sum, na.rm=T, fill=NA, align = 'left') / dt * 60
      sn_callrates[i,] <- rollapply(sn_timeLine[i,], width = dt, FUN = sum, na.rm=T, fill=NA, align = 'left') / dt * 60
      cc_prev_callrates[i,] <- rollapply(cc_timeLine[i,], width = dt_forage, FUN = sum, na.rm=T, fill = NA, align = 'right') / dt_forage * 60
      sn_prev_callrates[i,] <- rollapply(sn_timeLine[i,], width = dt_forage, FUN = sum, na.rm=T, fill = NA, align = 'right') / dt_forage * 60
      al_prev_callrates[i,] <- rollapply(al_timeLine[i,], width = dt_forage, FUN = sum, na.rm=T, fill = NA, align = 'right') / dt_forage * 60
    }
    
    #remove data for anything when audio was not on
    cc_callrates[!audioOn] <- NA
    sn_callrates[!audioOn] <- NA
    cc_prev_callrates[!audioOn] <- NA
    sn_prev_callrates[!audioOn] <- NA
    al_prev_callrates[!audioOn] <- NA
    
    #get number of callers for ccs and short notes in each time step
    cc_callers <- colSums(cc_prev_callrates > 0, na.rm=T)
    al_callers <- colSums(al_prev_callrates > 0, na.rm=T)
    tracked_callers <- colSums(!is.na(cc_prev_callrates)) #number of callers for which we have data at that timestep
    
    #fraction of individuals in the group giving cc's in the past dt_forage time
    frac_cc_callers <- cc_callers / tracked_callers
    frac_cc_callers[which(tracked_callers==0)] <- NA
    
    #Remove data when less than half the group gave a close call in the past dt_forage (60 sec) and when there was any alarm calling in the past dt_forage sec
    idxs_remove <- which(frac_cc_callers < 0.5 | al_callers > 0)
    cc_callrates[,idxs_remove] <- sn_callrates[,idxs_remove] <- NA
    
    #normalize call probability by individual - use z score
    cc_callrates_zscore <- sn_callrates_zscore <- matrix(NA, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      
      meancc <- mean(cc_callrates[i,], na.rm=T)
      sdcc <- sd(cc_callrates[i,],na.rm=T)
      cc_callrates_zscore[i,] <- (cc_callrates[i,] - meancc) / sdcc
      
      meansn <- mean(sn_callrates[i,], na.rm=T)
      sdsn <- sd(sn_callrates[i,],na.rm=T)
      sn_callrates_zscore[i,] <- (sn_callrates[i,] - meansn) / sdsn
      
    }
    
    #Get number of calls within radius R_neighbor in the past dt_forage time
    local_ccs <- local_sns <- matrix(NA, nrow = n_inds, ncol = n_times)
    for(i in 1:n_inds){
      neighbors_i <- dyad_dists[i,,] < R_neighbor
      local_ccs[i,] <- colSums(neighbors_i * cc_prev_callrates, na.rm=T)
      local_sns[i,] <- colSums(neighbors_i * sn_prev_callrates, na.rm=T)
      local_ccs[i,which(is.na(n_neighbors[i,]))] <- NA
      local_sns[i,which(is.na(n_neighbors[i,]))] <- NA
    }
    
    #STORE DATA
    data[[g]]$speeds <- speeds
    data[[g]]$cc_callrates <- cc_callrates
    data[[g]]$sn_callrates <- sn_callrates
    data[[g]]$cc_callrates_zscore <- cc_callrates_zscore
    data[[g]]$sn_callrates_zscore <- sn_callrates_zscore
    data[[g]]$n_neighbors <- n_neighbors
    data[[g]]$neighbor_change <- neighbor_change
    data[[g]]$n_tracked <- n_tracked
    data[[g]]$cc_prev_callrates <- cc_prev_callrates
    data[[g]]$al_prev_callrates <- al_prev_callrates
    data[[g]]$cc_callers <- cc_callers
    data[[g]]$al_callers <- al_callers
    data[[g]]$tracked_callers <- tracked_callers
    data[[g]]$surroundedness <- surroundedness
    data[[g]]$local_ccs <- local_ccs
    data[[g]]$local_sns <- local_sns
    
  }
}

#----------------------------------------------------------------------
#------------PLOTS------------
#----------------------------------------------------------------------

if(make_plots){
  for(g in 1:length(groupyears)){
    
    groupyear <- groupyears[[g]]
    
    speeds <- data[[g]]$speeds 
    cc_callrates <- data[[g]]$cc_callrates 
    sn_callrates <- data[[g]]$sn_callrates 
    cc_callrates_zscore <- data[[g]]$cc_callrates_zscore 
    sn_callrates_zscore <- data[[g]]$sn_callrates_zscore 
    neighbor_change <- data[[g]]$neighbor_change
    n_neighbors <- data[[g]]$n_neighbors
    n_tracked <- data[[g]]$n_tracked
    cc_callers <- data[[g]]$cc_callers
    sn_callers <- data[[g]]$sn_callers
    tracked_callers <- data[[g]]$tracked_callers
    surroundedness <- data[[g]]$surroundedness
    local_ccs <- data[[g]]$local_ccs
    local_sns <- data[[g]]$local_sns
    
    #bins
    n_neighbor_bins <- seq(-.5, max(n_neighbors,na.rm=T)+.5, 1)
    
    #----------Bivariate plots----
    #ROW 1: Individual speed vs mean call rates
    quartz(height = 14, width = 10)
    par(mfrow=c(4,2))
    #par(cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    mean_cc_rates <- mean_sn_rates <- mean_cc_rates_zscore <- mean_sn_rates_zscore <- rep(NA, length(speed_bins)-1)
    for(i in 1:(length(speed_bins)-1)){
      idxs <- which(speeds >= speed_bins[i] & speeds < speed_bins[i+1])
      mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
      mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
      mean_cc_rates_zscore[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
      mean_sn_rates_zscore[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
    }
    
    #make plots
    mids <- (speed_bins[1:(length(speed_bins)-1)] + speed_bins[2:length(speed_bins)]) / 2
    #plot(mids, mean_cc_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean CC rate (calls / min)', log = 'x')
    #plot(mids, mean_sn_rates, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean CC rate (calls / min)', log = 'x')
    plot(mids, mean_cc_rates_zscore, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean CC z-score', log = 'x')
    abline(h=0, lty = 2)
    plot(mids, mean_sn_rates_zscore, pch = 19, xlab = 'Speed (m / min)', ylab = 'Mean SN z-score', log = 'x')
    abline(h=0, lty = 2)
    
    #ROW 2: N neighbors vs mean call rates
    mean_cc_rates <- mean_sn_rates <- mean_cc_rates_zscore <- mean_sn_rates_zscore <- rep(NA, length(n_neighbor_bins)-1)
    for(i in 1:(length(n_neighbor_bins)-1)){
      idxs <- which(n_neighbors >= n_neighbor_bins[i] & n_neighbors < n_neighbor_bins[i+1])
      mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
      mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
      mean_cc_rates_zscore[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
      mean_sn_rates_zscore[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
    }
    
    #make plot
    mids <- (n_neighbor_bins[1:(length(n_neighbor_bins)-1)] + n_neighbor_bins[2:length(n_neighbor_bins)])/2
    #plot(mids, mean_cc_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean CC rate (calls / min)')
    #plot(mids, mean_sn_rates, pch = 19, xlab = '# neighbors', ylab = 'Mean SN rate (calls / min)')
    plot(mids, mean_cc_rates_zscore, pch = 19, xlab = '# neighbors', ylab = 'Mean CC z-score')
    abline(h=0, lty = 2)
    plot(mids, mean_sn_rates_zscore, pch = 19, xlab = '# neighbors', ylab = 'Mean SN z-score')
    abline(h=0, lty = 2)
    
    #ROW 3: surroundedness vs call rate
    mean_cc_rates <- mean_sn_rates <- mean_cc_rates_zscore <- mean_sn_rates_zscore <- rep(NA, length(surroundedness_bins)-1)
    for(i in 1:(length(surroundedness_bins)-1)){
      idxs <- which(surroundedness >= surroundedness_bins[i] & surroundedness < surroundedness_bins[i+1])
      mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
      mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
      mean_cc_rates_zscore[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
      mean_sn_rates_zscore[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
    }
    
    #make plot
    mids <- (surroundedness_bins[1:(length(surroundedness_bins)-1)] + surroundedness_bins[2:length(surroundedness_bins)])/2
    #plot(mids, mean_cc_rates, pch = 19, xlab = 'Surroundedness', ylab = 'Mean CC rate')
    #plot(mids, mean_sn_rates, pch = 19, xlab = 'Surroundedness', ylab = 'Mean SN rate')
    plot(mids, mean_cc_rates_zscore, pch = 19, xlab = 'Surroundedness', ylab = 'Mean CC z-score')
    abline(h=0, lty = 2)
    plot(mids, mean_sn_rates_zscore, pch = 19, xlab = 'Surroundedness', ylab = 'Mean SN z-score')
    abline(h=0, lty = 2)
    
    #ROW 4: local calls in past dt_forage vs call rate
    mean_cc_rates <- mean_sn_rates <- mean_cc_rates_zscore <- mean_sn_rates_zscore <- rep(NA, length(local_call_bins)-1)
    for(i in 1:(length(local_call_bins)-1)){
      idxs <- which(local_ccs >= local_call_bins[i] & local_ccs < local_call_bins[i+1])
      mean_cc_rates[i] <- mean(cc_callrates[idxs], na.rm=T)
      mean_cc_rates_zscore[i] <- mean(cc_callrates_zscore[idxs], na.rm=T)
      idxs <- which(local_sns >= local_call_bins[i] & local_sns < local_call_bins[i+1])
      mean_sn_rates[i] <- mean(sn_callrates[idxs], na.rm=T)
      mean_sn_rates_zscore[i] <- mean(sn_callrates_zscore[idxs], na.rm=T)
    }
    
    #make plot
    mids <- (local_call_bins[1:(length(local_call_bins)-1)] + local_call_bins[2:length(local_call_bins)])/2
    #plot(mids, mean_cc_rates, pch = 19, xlab = 'Local call rate', ylab = 'Mean CC rate', log = 'x')
    #plot(mids, mean_sn_rates, pch = 19, xlab = 'Local call rate', ylab = 'Mean SN rate', log = 'x')
    plot(mids, mean_cc_rates_zscore, pch = 19, xlab = 'Local call rate', ylab = 'Mean CC z-score', log = 'x')
    abline(h=0, lty = 2)
    plot(mids, mean_sn_rates_zscore, pch = 19, xlab = 'Local call rate', ylab = 'Mean SN z-score', log = 'x')
    abline(h=0, lty = 2)
    
    dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_neighbors_zscore_',groupyear,'.pdf'))
    dev.off()
    
    #---------------3d plots--------------------------
    quartz(width = 8, height = 14)
    par(mfrow=c(4,2))
    #speed and number of neighbors
    cc_rates_2d <- sn_rates_2d <- cc_rates_2d_zscore <- sn_rates_2d_zscore <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(n_neighbor_bins)-1)
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
          cc_rates_2d_zscore[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
          sn_rates_2d_zscore[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
        }
      }
    }
    image.plot(z=cc_rates_2d_zscore, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main = 'CC')
    image.plot(z=sn_rates_2d_zscore, col = viridis(256), xlab = 'Speed ', ylab ='# neighbors',main='SN')
    
    #surroundedness and speed
    cc_rates_2d <- sn_rates_2d <- cc_rates_2d_zscore <- sn_rates_2d_zscore <- nsamps <- matrix(NA, nrow = length(speed_bins)-1, ncol = length(surroundedness_bins)-1)
    for(i in 1:(length(speed_bins)-1)){
      for(j in 1:(length(surroundedness_bins)-1)){
        idxs <- which(speeds >= speed_bins[i] &
                        speeds < speed_bins[i+1] &
                        surroundedness >= surroundedness_bins[j] &
                        surroundedness < surroundedness_bins[j+1])
        nsamps[i,j] <- length(idxs)
        if(length(idxs)>1000){
          cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
          sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
          cc_rates_2d_zscore[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
          sn_rates_2d_zscore[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
        }
      }
    }
    
    image.plot(z=cc_rates_2d_zscore, col = viridis(256), xlab = 'Speed ', ylab ='Surroundedness',main = 'CC')
    image.plot(z=sn_rates_2d_zscore, col = viridis(256), xlab = 'Speed ', ylab ='Surroundedness',main='SN')
    
    #surroundedness and number of neighbors
    cc_rates_2d <- sn_rates_2d <- cc_rates_2d_zscore <- sn_rates_2d_zscore <- nsamps <- matrix(NA, nrow = length(surroundedness_bins)-1, ncol = length(n_neighbor_bins)-1)
    for(i in 1:(length(surroundedness_bins)-1)){
      for(j in 1:(length(n_neighbor_bins)-1)){
        idxs <- which(surroundedness >= surroundedness_bins[i] &
                        surroundedness < surroundedness_bins[i+1] &
                        n_neighbors >= n_neighbor_bins[j] &
                        n_neighbors < n_neighbor_bins[j+1])
        nsamps[i,j] <- length(idxs)
        if(length(idxs)>1000){
          cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
          sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
          cc_rates_2d_zscore[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
          sn_rates_2d_zscore[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
        }
      }
    }
    
    image.plot(z=cc_rates_2d_zscore, col = viridis(256), xlab = 'Surroundedness ', ylab ='# neighbors',main = 'CC')
    image.plot(z=sn_rates_2d_zscore, col = viridis(256), xlab = 'Surroundedness ', ylab ='# neighbors',main='SN')
    
    #----------Call rate vs number of neighbors and their call rates----------
    #surroundedness and number of neighbors
    cc_rates_2d <- sn_rates_2d <- cc_rates_2d_zscore <- sn_rates_2d_zscore <- nsamps <- matrix(NA, nrow = length(surroundedness_bins)-1, ncol = length(n_neighbor_bins)-1)
    for(i in 1:(length(local_call_bins)-1)){
      for(j in 1:(length(n_neighbor_bins)-1)){
        
        #CC
        idxs <- which(local_ccs >= local_call_bins[i] &
                        local_ccs < local_call_bins[i+1] &
                        n_neighbors >= n_neighbor_bins[j] &
                        n_neighbors < n_neighbor_bins[j+1])
        nsamps[i,j] <- length(idxs)
        if(length(idxs)>1000){
          cc_rates_2d[i,j] <- mean(cc_callrates[idxs],na.rm=T)
          cc_rates_2d_zscore[i,j] <- mean(cc_callrates_zscore[idxs],na.rm=T)
        }
        
        #SN
        idxs <- which(local_sns >= local_call_bins[i] &
                        local_sns < local_call_bins[i+1] &
                        n_neighbors >= n_neighbor_bins[j] &
                        n_neighbors < n_neighbor_bins[j+1])
        
        if(length(idxs)>1000){
          sn_rates_2d[i,j] <- mean(sn_callrates[idxs],na.rm=T)
          sn_rates_2d_zscore[i,j] <- mean(sn_callrates_zscore[idxs],na.rm=T)
        }
        
        
      }
    }
    
    image.plot(z=cc_rates_2d_zscore, col = viridis(256), xlab = 'Local call rate ', ylab ='# neighbors',main = 'CC')
    image.plot(z=sn_rates_2d_zscore, col = viridis(256), xlab = 'Local call rate ', ylab ='# neighbors',main='SN')
    
    dev.copy2pdf(file = paste0(outdir, 'callrate_vs_speed_neighbors_surroundedness_callrate_3d_',groupyear,'.pdf'))
    dev.off()
    
    
  }
}

#Look at speed results when at least half group giving close calls (filter out non-foraging periods) - apply to all analyses (50% of group within past 60 sec, no alarms in past 60 sec) - done
#Call rate vs. number of neighbors and surroundedness (circular variance) - done
#Call rate vs number of neighbors and their call rates (?)
#Call-response-response triplets: AAA AAB ABA ABC ABB
#ERROR BARS - bootstrap by day
