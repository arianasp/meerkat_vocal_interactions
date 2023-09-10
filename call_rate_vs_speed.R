#CALL RATE VS SPEED

#This script computes the mean rate of calling as a function of individual speed 
#and number of neighbors, and produces plots.
#The script also has some code to compute other metrics, however these are not 
#addressed in the current manuscript (maybe will be in future papers).

#Libraries
library(lubridate)
library(zoo)
library(fields)
library(viridis)
library(lme4)

#-------------FLAGS-------------
compute <- F #whether to perform computations and save (if not, they will be loaded from a precomputed file)
make_plots <- T #whether to produce plots or not

#--------------PARAMS-------------
#group years
groupyears <- c('HM2017','HM2019','L2019') 

#time window (future) over which to compute call rate and speed
dt <- 10 

#input and output directories
indir <- '/Users/Ari/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/'
outdir <- '/Users/Ari/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/plots_2023-08-30_afterccfilt/'

#name of file (csv) where call data is stored (inside indir)
audio.file <- 'all_calls_sync_resolved_2023-03-23_cc_sn_filt.csv'

#minimum number of individuals tracked in order to include data
min_n_tracked <- 5

#local neighborhood radius (in meters) - for computing number of neighbors
R_neighbor <- 5 

#time window to use to specfiy whether the group is foraging (if at least half of the group gives a cc during this window), and other metrics related to the past behavior of the group
dt_past <- 60 

#call types to use (they will be grepl-ed for)
#for paper, need cc, sn and al (al just for filtering out data when alarms occurred)
calltypes <- c('cc','sn','al') 

#minimum data per individual to plot a curve for them (only used in plotting code)
min_dat <- 3600*2 #at least 2 hrs

#-----------MAIN-----------

#-------PERFORM COMPUTATIONS OF METRICS IF REQRUIED--------
if(compute){
  
  #STORE PARAMETERS IN PARAMS OBJECT
  params <- list()
  params$groupyears <- groupyears
  params$dt <- dt
  params$indir <- indir
  params$outdir <- outdir
  params$min_n_tracked <- min_n_tracked
  params$R_neighbor <- R_neighbor
  params$dt_past <- dt_past
  params$calltypes <- calltypes

  #------LOOP OVER GROUPS AND STORE METRICS DATA------
  data <- list()
  
  setwd(indir)
  for(g in 1:length(groupyears)){
    
    groupyear <- groupyears[g]
    data[[g]] <- list()
    
    print(groupyear)
    
    #--------LOAD DATA----------
    load(paste0(groupyear,'_COORDINATES_all_sessions.RData'))
    calls <- read.csv(file = paste0(audio.file), header=T, stringsAsFactors =F, sep = ',')
    load(paste0(groupyear,'_DATAPRESENCE_all_sessions.RData'))
    
    #get data associated with that groupyear - store with general names
    timeline <- eval(as.name(paste(groupyear,'timeLine', sep = '_')))
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
    
    #get calls only for that groupyear
    dates <- unique(as.Date(timeline))
    dates <- gsub('-','',dates)
    calls <- calls[which(calls$date %in% dates),]
    
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
    n_neighbors[,which(n_tracked < min_n_tracked)] <- NA #remove data for all if not enough inds tracked
    n_neighbors[!gpsOn] <- NA #remove data for any individual who is not currently tracked
    
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
    
    #get datetime and time index for each call in the calls table
    calls$datetime <- as.POSIXct(calls$t0GPS_UTC, tz = 'UTC')
    calls$tidx <- match(as.character(calls$datetime),as.character(as.POSIXct(timeline,tz='UTC')))
    
    #calculate number of calls of each type in each second, for each individual
    call_timelines <- list()
    for(c in 1:length(calltypes)){
      call_timelines[[calltypes[c]]] <- matrix(0, nrow = n_inds, ncol = n_times)
    }
    
    for(c in 1:length(calltypes)){
      for(i in 1:n_inds){
        call_times <- calls$tidx[which(grepl(calltypes[c],calls$stn_call_type) 
                                      #& calls$pred_focalType=='F' #all calls are now focal in the new dataset
                                      & calls$ind == indInfo$code[i])]
        call_times <- call_times[!is.na(call_times)]
        for(t in 1:length(call_times)){
          call_timelines[[calltypes[c]]][i,call_times[t]] <- call_timelines[[calltypes[c]]][i,call_times[t]] + 1
        }
      }
      #remove data when audio was not labeled
      call_timelines[[calltypes[c]]][!audioOn] <- NA
    }
    
    #get call rates in the time window dt (future) and callrates in the past dt_past
    callrates <- callrates_past <- list()
    for(c in 1:length(calltypes)){
      callrates[[calltypes[c]]] <- callrates_past[[calltypes[c]]] <- matrix(0, nrow = n_inds, ncol = n_times)
      for(i in 1:n_inds){
        callrates[[calltypes[c]]][i,] <- rollapply(call_timelines[[calltypes[c]]][i,], width = dt, FUN = sum, na.rm=T, fill=NA, align = 'left') / dt * 60
        callrates_past[[calltypes[c]]][i,] <- rollapply(call_timelines[[calltypes[c]]][i,], width = dt_past, FUN = sum, na.rm=T, fill=NA, align = 'right') / dt_past * 60
      }
      #remove data for anything when audio was not on
      callrates[[calltypes[c]]][!audioOn] <- NA
      callrates_past[[calltypes[c]]][!audioOn] <- NA
    }
    
    #Removing data from times when group was not foraging OR was alarming, in the past dt_past time
    #get number of cc callers
    cc_callers <- colSums(callrates_past[['cc']] > 0, na.rm=T)
    tracked_callers <- colSums(!is.na(callrates_past[['cc']]))
    frac_cc_callers <- cc_callers / tracked_callers
    frac_cc_callers[which(tracked_callers < min_n_tracked)] <- NA
    
    #number of alarm calls in the past dt_forage
    alarm_rate <- colSums(callrates_past[['al']], na.rm=T)
    alarm_rate[which(tracked_callers == 0)] <- NA
    
    #Remove data when less than half the group gave a close call in the past dt_forage (60 sec) and when there was any alarm calling in the past dt_forage sec
    idxs_remove <- which(frac_cc_callers < 0.5 | alarm_rate > 0)
    for(c in 1:length(calltypes)){
      callrates[[calltypes[c]]][,idxs_remove] <- NA
    }
    
    #normalize call probability by individual - use z score
    callrates_zscore <- callrates
    for(c in 1:length(calltypes)){
      for(i in 1:n_inds){
        meanrate <- mean(callrates[[calltypes[c]]][i,], na.rm=T)
        sdrate <- sd(callrates[[calltypes[c]]][i,], na.rm=T)
        callrates_zscore[[calltypes[c]]][i,] <- (callrates[[calltypes[c]]][i,] - meanrate) / sdrate
      }
    }
    
    #Get number of calls within radius R_neighbor in the past dt_past time
    callrates_local <- callrates
    for(c in 1:length(calltypes)){
      for(i in 1:n_inds){
        neighbors_i <- dyad_dists[i,,] < R_neighbor
        callrates_local[[calltypes[c]]][i,] <- colSums(neighbors_i * callrates_past[[calltypes[c]]], na.rm=T)
        callrates_local[[calltypes[c]]][i,which(is.na(n_neighbors[i,]))] <- NA
      }
    }
    
    #STORE DATA
    data[[g]]$calltypes <- calltypes
    data[[g]]$n_tracked_gps <- n_tracked
    data[[g]]$n_tracked_audio <- tracked_callers
    data[[g]]$speeds <- speeds
    data[[g]]$n_neighbors <- n_neighbors
    data[[g]]$neighbor_change <- neighbor_change
    data[[g]]$surroundedness <- surroundedness
    data[[g]]$callrates_fut <- callrates
    data[[g]]$callrates_past <- callrates_past
    data[[g]]$callrates_zscore <- callrates_zscore
    data[[g]]$callrates_local <- callrates_local
    data[[g]]$timeline <- timeline 
  }
  
  save(list = c('params','data'), file = paste0(params$outdir,'metrics_data.RData'))
  
}

#---------CREATE PLOTS----------
if(make_plots){
  
  #Load data
  setwd(outdir)
  load(paste0(outdir,'metrics_data.RData'))
  
  #aggregate all data into single vectors for aggregated measures (across all groups and individuals)
  speeds_all <- nbrs_all <- cc_rates_zscore_all <- sn_rates_zscore_all <- boot_chunks_all <- c()
  for(g in 1:length(groupyears)){
    speeds_all <- c(speeds_all, c(data[[g]]$speeds))
    nbrs_all <- c(nbrs_all, c(data[[g]]$n_neighbors))
    cc_rates_zscore_all <- c(cc_rates_zscore_all, c(data[[g]]$callrates_zscore[['cc']]))
    sn_rates_zscore_all <- c(sn_rates_zscore_all, c(data[[g]]$callrates_zscore[['sn']]))
    boot_chunks_all <- c(boot_chunks_all, as.Date(c(data[[g]]$timeline)))
  }
  
  #get total number of individuals across all groupyears
  n_inds_tot <- 0
  for(g in 1:length(groupyears)){
    n_inds_tot <- n_inds_tot + nrow(data[[g]]$speeds)
  }
  
  #get speed bins (10% quantiles, 10 bins)
  speed_bins <- quantile(speeds_all, seq(0,1,length.out=11), na.rm=T)
  
  #get neighbor bins (0, 1, 2, 3, >=4)
  n_neighbor_bins <- c(seq(0,4,1), max(nbrs_all,na.rm=T))
  
  #------INITIALIZE PLOTTING-----
  quartz(height = 7, width = 7)
  par(mfcol = c(2,2))
  par(cex.lab=1.5)
  par(mar = c(5,5,1,1))
  
  #-----COLLECT DATA FOR PLOTTING---
  
  #initialize empty matrices for individual-level (row = individual, col = speed or nbr bin)
  cc_vs_speed_ind <- sn_vs_speed_ind <- matrix(NA, nrow = n_inds_tot, ncol = length(speed_bins)-1)
  cc_vs_nbr_ind <- sn_vs_nbr_ind <- matrix(NA, nrow = n_inds_tot, ncol = length(n_neighbor_bins)-1)
  
  #initiative empty vector for aggregate (across all inds and groupyears) lines
  cc_vs_speed_all <- sn_vs_speed_all <- rep(NA, length(speed_bins)-1)
  cc_vs_nbr_all <- sn_vs_nbr_all <- rep(NA, length(n_neighbor_bins)-1)
  
  #----Call rate vs speed plots----
  bins <- speed_bins #bins to use 
  mids <- (bins[1:(length(bins)-1)] + bins[2:length(bins)])/2 #midpoints of the bins (for plotting)
  
  #loop over groupyears and collect data into matrices and vectors
  start_idx <- 0 #index of where to start storing data for each group in the matrices
  for(g in 1:length(groupyears)){ #for each group...
    n_inds <- nrow(data[[g]]$speeds) #get number of individuals in the group
    idxs <- start_idx + 1:n_inds #get indexes in the matrices where data will be stored
    for(i in 1:n_inds){ #for each individual...
      cc_curr <- data[[g]]$callrates_zscore[['cc']][i,] #get cc data for that ind
      sn_curr <- data[[g]]$callrates_zscore[['sn']][i,] #get sn data for that ind
      xcurr <- data[[g]]$speeds[i,] #get the data for the independent variable (here speed)
      if(sum(!is.na(xcurr) &!is.na(cc_curr))>min_dat){ #include only if there is enoough data for that ind (must have both gps and audio data)
        for(j in 1:(length(bins)-1)){ #loop over bins...
          idxs_curr <- which(xcurr >= bins[j] & xcurr < bins[j+1]) #get indexes where xcurr (speed) fell within a given bin
          cc_vs_speed_ind[idxs[i], j] <- mean(cc_curr[idxs_curr], na.rm=T) #compute mean cc rate for those indexes
          sn_vs_speed_ind[idxs[i], j] <- mean(sn_curr[idxs_curr], na.rm=T) #compute mean sn rate for those indexes
        }
      }
    }
    start_idx <- start_idx + n_inds #update start index to store the next group in the next set of rows in the output matrix
  }
  
  #compute means at the group level and store in the output vector
  for(j in 1:(length(bins)-1)){ #for each bins
    curr_idxs <- which(speeds_all >= bins[j] & speeds_all < bins[j+1]) #get indexes associated with speeds in that bin
    cc_vs_speed_all[j] <- mean(cc_rates_zscore_all[curr_idxs],na.rm=T) #compute mean cc rate
    sn_vs_speed_all[j] <- mean(sn_rates_zscore_all[curr_idxs],na.rm=T) #compute mean sn rate
  }
  
  #make plots for call rate vs speed
  #cc plot
  plot(NULL, xlim = c(min(mids), max(mids)), ylim = c(-.5,.5), log = 'x', xlab = 'Speed (m / min)', ylab = 'Call rate (z-score)')
  for(i in 1:n_inds_tot){
    lines(mids, cc_vs_speed_ind[i,], col = '#FF000044')
  }
  lines(mids, cc_vs_speed_all, lwd = 3, col = 'red')
  points(mids, cc_vs_speed_all, cex = 1.5, pch = 19, col = 'red')
  legend('topleft','Close calls', bty='n', cex = 1.2, text.col = 'red')
  abline(h=0,lty=2)
  
  #sn plot
  plot(NULL, xlim = c(min(mids), max(mids)), ylim = c(-.5,.5), log = 'x', xlab = 'Speed (m / min)', ylab = 'Call rate (z-score)')
  for(i in 1:n_inds_tot){
    lines(mids, sn_vs_speed_ind[i,], col = '#0000FF44')
  }
  lines(mids, sn_vs_speed_all, lwd = 3, col = 'blue')
  points(mids, sn_vs_speed_all, cex = 1.5, pch = 19, col = 'blue')
  legend('topleft','Short note calls', bty='n', cex = 1.2, text.col = 'blue')
  abline(h=0,lty=2)
  
  #-----Call rate vs n nbrs plots----
  #same as above, but for a different independent variable (number of neighbors)
  bins <- n_neighbor_bins #get appropriate bins
  start_idx <- 0
  
  #loop over groupyears and collect data
  for(g in 1:length(groupyears)){
    n_inds <- nrow(data[[g]]$n_neighbors)
    idxs <- start_idx + 1:n_inds
    for(i in 1:n_inds){
      cc_curr <- data[[g]]$callrates_zscore[['cc']][i,]
      sn_curr <- data[[g]]$callrates_zscore[['sn']][i,]
      xcurr <- data[[g]]$n_neighbors[i,] #use number of neighbors as independent variable
      if(sum(!is.na(xcurr) &!is.na(cc_curr))>min_dat){
        for(j in 1:(length(bins)-1)){
          idxs_curr <- which(xcurr >= bins[j] & xcurr < bins[j+1])
          cc_vs_nbr_ind[idxs[i], j] <- mean(cc_curr[idxs_curr], na.rm=T)
          sn_vs_nbr_ind[idxs[i], j] <- mean(sn_curr[idxs_curr], na.rm=T)
        }
      }
    }
    start_idx <- start_idx + n_inds
  }
  
  #same at group level
  for(j in 1:(length(bins)-1)){
    curr_idxs <- which(nbrs_all >= bins[j] & nbrs_all < bins[j+1])
    cc_vs_nbr_all[j] <- mean(cc_rates_zscore_all[curr_idxs],na.rm=T)
    sn_vs_nbr_all[j] <- mean(sn_rates_zscore_all[curr_idxs],na.rm=T)
  }
  
  #plot
  #cc plot
  plot(NULL, xlim = c(0, bins[length(bins)-1]), ylim = c(-.5,.5), xlab = 'Number of neighbours', ylab = 'Call rate (z-score)', xaxt= 'n')
  axis(1, at = 0:4, labels = c('0','1','2','3','4+'))
  for(i in 1:n_inds_tot){
    lines(bins[1:(length(bins)-1)], cc_vs_nbr_ind[i,], col = '#FF000044')
  }
  lines(bins[1:(length(bins)-1)], cc_vs_nbr_all, lwd = 3, col = 'red')
  points(bins[1:(length(bins)-1)], cc_vs_nbr_all, cex = 1.5, pch = 19, col = 'red')
  legend('topleft','Close calls', bty='n', cex = 1.2, text.col = 'red')
  abline(h=0,lty=2)
  
  #sn plot
  plot(NULL, xlim = c(0, bins[length(bins)-1]), ylim = c(-.5,.5), xlab = 'Number of neighbours', ylab = 'Call rate (z-score)', xaxt='n')
  axis(1, at = 0:4, labels = c('0','1','2','3','4+'))
  for(i in 1:n_inds_tot){
    lines(bins[1:(length(bins)-1)], sn_vs_nbr_ind[i,], col = '#0000FF44')
  }
  lines(bins[1:(length(bins)-1)], sn_vs_nbr_all, lwd = 3, col = 'blue')
  points(bins[1:(length(bins)-1)], sn_vs_nbr_all, cex = 1.5, pch = 19, col = 'blue')
  legend('topleft','Short note calls', bty='n', cex = 1.2, text.col = 'blue')
  abline(h=0,lty=2)
  
  dev.copy2pdf(file = paste0(outdir,'callrate_vs_speed_and_neighbors.pdf'))
  
}
