#Spatiotemporal call clustering - do calls cluster together in space and time?

#This code runs the call spatiotemporal clustering analysis described in the paper 
#"Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups"
#It starts from the raw GPS and call label data and outputs processed data that are used in a later plotting 
#script (plot_spatiotemporal_call_clustering.R)

#How does it work?
#First, we define a clustering metric.
#Here we use modification of the "Knox statistic" (actually kind of numerical derivative of the Knox statistic, described below)
#As a metric of clustering, we compute the number of pairs of calls that 
#fall within a distance range r to r + dr and a temporal range t to t + dt. We then noramlize this value by the total
#number of pairs of calls within all distance and time ranges

#Next, we create a null model
#Because meerkats move in cohesive groups, there will trivially be some spatial and temporal clustering of calls.
#Just because the meerkats themselves are spatially and temporally clustered. So we want to see whether the call
#distribution is MORE clustered than the meerkat distribution. To test this, we break the link between meerkat
#trajectories and their calling behavior. We generate "pseudo-meerkats" which have the calling behavior of one 
#individual but the movement behavior of another by shuffling the identities associated with the gps trajectories
#within each day. We then recompute the clustering metric described above.

#Finally, we compare data to null model
#To do so, we take the ratio of the clustering statistic K in the real vs permuted data. We then plot the log of this ratio.

#----------LIBRARIES--------
library(stringr)
library(lubridate)

#--------------PARAMETERS----------

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------
#directories where data is stored (these are probably the same directory)
#audiodir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 
#gpsdir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 
audiodir <- '~/Desktop/meerkat_data_anon/'
gpsdir <- '~/Desktop/meerkat_data_anon/'

#directory where you would like to store the processed data (can be the same directory too if you like)
#savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/precomputed/'
savedir <- '~/Desktop/meerkat_data_anon/precomputed/'

#call type
callType <- 'sn'

#test on a smaller number of bins and with only one null model?
testflag <- F

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------

#time windows (sec)
time.windows <- c(0,1,2,5,10,20,30,60,120,300,600,1800,3600)

#distance windows (m)
dist.windows <- c(0,2,5,10,25,50,100,200)

#number of randomizations
n.rands <- 100

#list of sessions to use
sessions <- c('HM2017','HM2019','L2019')

#minimum number of individual present to include in dataset
min.inds.present <- 5

#if testing, modify parameters
if(testflag){
  time.windows <- c(1,600,3600)
  dist.windows <- c(2,5,50)
  n.rands <- 1
}

#-------- FUNCTIONS ------

#Compute the modified Knox index (see description above)
#INPUTS:
# calls.include: call data frame
# allX, allY: matrices [n.inds x n.times] of x and y positions respectively
# dist.windows: distance windows used by the metric (m)
# time.windows: time windows used by the metric (sec)
# diff.inds.only: boolean variable determining whether we should remove same-individual pairs of calls (T) or not (F)
#OUTPUTS:
# out: a list containing
#   out$dist.windows: distance window (same as input)
#   out$time.windows: time window (same as input)
#   out$K: matrix of the modified Knox K metric (see explanation above)
#   out$npairs.dist.time: matrix of number of pairs of calls within a given distance range and time range (numerator of K)
#   out$tot.pairs: matrix of total pairs of calls in all distance and time ranges (denominator of K)
compute_Knox_indexes <- function(calls.include, allX, allY, dist.windows, time.windows, diff.inds.only = F){
  
  #get x and y positions at the times of calls
  x_call <- allX[cbind(calls.include$ind.idx, calls.include$time.idx)]
  y_call <- allY[cbind(calls.include$ind.idx, calls.include$time.idx)]
  t_call <- as.numeric(as.POSIXlt(calls.include$t0GPS_UTC))
  
  #remove NAs
  non.na.idxs <- which(!is.na(x_call) & !is.na(y_call) & !is.na(t_call))
  if(length(non.na.idxs)==0){
    stop('All values are NAs')
  }
  x_call <- x_call[non.na.idxs]
  y_call <- y_call[non.na.idxs]
  t_call <- t_call[non.na.idxs]
  
  #construct a spatial distance matrix (actually formatted as a vector) for the dist between every pair of calls
  spatial.dist <- dist(cbind(x_call, y_call))
  
  #construct a temporal distance matrix for the temporal distance between every pair of calls
  temporal.dist <- dist(t_call)
  
  #remove pairs of calls given by the same individual if specified
  if(diff.inds.only){
    #construct matrix to determine if the caller was the same caller (1) or a different (0)
    diff.inds <- dist(calls.include$ind.idx) != 0
    
    #remove same caller data
    spatial.dist <- spatial.dist[diff.inds]
    temporal.dist <- temporal.dist[diff.inds]
  }
  
  #matrices to store output
  npairs.dist.time <- matrix(NA, nrow = length(dist.windows)-1, ncol = length(time.windows)-1)
  for(r in 2:length(dist.windows)){
    for(t in 2:length(time.windows)){
      
      #get booleans indicating whether a given call pair falls within a given distance bin and time bin
      in.dist.bin <- (spatial.dist >= dist.windows[r-1]) & (spatial.dist < dist.windows[r])
      in.time.bin <- (temporal.dist >= time.windows[t-1]) & (temporal.dist < time.windows[t])
      
      #get the total number of calls within a distance window and a time window
      npairs.dist.time[r-1,t-1] <- sum(in.dist.bin & in.time.bin, na.rm=T)
    }
  }
  
  #normalize by total number of pairs of calls within all time and distance bins
  tot.pairs <- sum(npairs.dist.time, na.rm=T)
  K <- npairs.dist.time / tot.pairs
  
  #list for output
  out <- list()
  out$dist.windows <- dist.windows
  out$time.windows <- time.windows
  out$K <- npairs.dist.time / tot.pairs
  out$npairs.dist.time <- npairs.dist.time
  out$tot.pairs <- tot.pairs
  
  return(out)
}

#function to simplify the name of certain objects in the environment (removing "pat" argument from their names) if reverse=F
#OR adding "pat" argument at the beginning of all object provided in argument "names" if reverse=T
#returns a vector of the new names
simplifyNames <- function(pat,reverse=F,items=ls(name=.GlobalEnv, pattern=pat)){
  
  for(item in items){
    if(reverse){
      assign(paste(pat,item,sep=""),get(item,envir=.GlobalEnv),envir=.GlobalEnv)
    }else{
      assign(str_remove(item,pat),get(item,envir=.GlobalEnv),envir=.GlobalEnv)
    }
    rm(list=item,envir=.GlobalEnv)
  }
  if(reverse){
    return(paste(pat,items,sep=""))
  }else{
    return(str_remove(items,pat))
  }
}

#------------SETUP-------------

print('Running call clustering analysis')

#load audio data
calls.allsessions <- read.csv(paste0(audiodir,'all_calls_sync_resolved_2023-09-10_cc_sn_filt_with_oor_anon.csv' ), header=T, stringsAsFactors=F)

timestamp()
for(sess.idx in 1:length(sessions)){

  #------------LOAD DATA -------
  session <- sessions[sess.idx]
  
  print(session)
  
  #loading spatial data and making names easier
  longNames <- c(load(paste0(gpsdir,session,"_COORDINATES_all_sessions_anon.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
  
  #get dates associated with that session
  dates <- unique(date(timeLine))
  dates <- gsub('-','',dates) #convert to same format as in audio label table
  
  #get individuals associated with that session
  inds <- indInfo$code
  
  #subset to calls in that session
  calls.all <- calls.allsessions[which(calls.allsessions$date %in% dates & calls.allsessions$ind %in% inds),]
  
  #add individual index
  calls.all$ind.idx <- match(calls.all$ind,indInfo$code)
  
  #add time index (note: rounds down the prior second, shouldn't matter for this analysis)
  call.times <- sub('\\..*','',calls.all$t0GPS_UTC)
  calls.all$time.idx <- match(call.times, timeLine)
  
  #get number of individuals
  n.inds <- nrow(indInfo)
  
  #get number of time steps
  n.times <- ncol(allX)
  
  #get beginning and end times of tracking periods (both GPS and audio) as well as which individuals have both GPS and audio
  #Note: This will not account for instances where there are gaps in audio labeling (but this should be rare, unlikely to systematically affect this analysis as it will be same in real vs randomized data)
  labeled.intervals <- data.frame(date = dates, 
                               audio.start = rep(NA,length(dates)),
                               audio.end = rep(NA, length(dates)))
  inds.present <- matrix(FALSE, nrow = nrow(labeled.intervals), ncol = n.inds)
  for(date.idx in 1:length(dates)){
    
    #get all calls on that date
    calls.date <- calls.all[which(calls.all$date == dates[date.idx]),]
    
    #get the individuals that have labeled audio on that date
    inds.date <- unique(calls.date$ind.idx)
    inds.present[date.idx, inds.date] <- TRUE
    start.times <- end.times <- rep(NA, length(inds.date))
    for(ind.idx in 1:length(inds.date)){
      ind <- inds.date[ind.idx]
      
      #find indexes to start and end marker(s) for that individual
      start.idxs <- which(calls.date$entryName %in% c('start','START') & calls.date$ind.idx == ind)
      end.idxs <- which(calls.date$entryName %in% c('end','END','stop','STOP') & calls.date$ind.idx == ind)
      
      #get times associated with them
      if(length(start.idxs)>0){
        start.times.ind <- calls.date$time.idx[start.idxs]
      } else{
        #if there were no 'start' labels found, use the earliest time of any entry for that individual
        #this gives warnings if there are no calls, but we deal with this later so suppress them
        start.times.ind <- suppressWarnings(min(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T))
      }
      if(length(end.idxs)>0){
        end.times.ind <- calls.date$time.idx[end.idxs]
      } else{
        #if there were no 'start' labels found, use the earliest time of any entry for that individual
        #this gives warnings if there are no calls, but we deal with this later so suppress them
        end.times.ind <- suppressWarnings(max(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T))
      }
      
      #if there is more than one start label, use the latest one
      start.time.ind <- max(start.times.ind)
      
      #if there is a more than one end label, use the earliest one
      end.time.ind <- min(end.times.ind)
      
      #if the start or end time index is NA, this means it extends outside the bounds of the GPS data - replace it with the first / last non-NA time index point
      #this gives warnings if there are no calls, but we deal with this later so suppress them
      start.time.ind <- suppressWarnings(min(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T))
      end.time.ind <- suppressWarnings(max(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T))
      
      #store the data for that individual's start and end time
      start.times[ind.idx] <- start.time.ind
      end.times[ind.idx] <- end.time.ind
    }
    
    #get the latest start time and earliest end time across all individuals - this is the recording interval
    latest.start <- max(start.times, na.rm=T)
    earliest.end <- min(end.times, na.rm=T)
    labeled.intervals$audio.start[date.idx] <- latest.start
    labeled.intervals$audio.end[date.idx] <- earliest.end
  }
  
  #if there are any intervals that don't contain any time, remove from labeled.intervals and inds.present, and output a message
  empty.intervals <- which(is.infinite(labeled.intervals$audio.start) | labeled.intervals$audio.start >= labeled.intervals$audio.end)
  if(length(empty.intervals)>0){
    print(paste('The interval for date', dates[empty.intervals], 'ends before it starts or did not contain audio data - removing'))
    labeled.intervals <- labeled.intervals[-empty.intervals,]
    inds.present <- inds.present[-empty.intervals,]
    dates <- dates[-empty.intervals]
  }
  
  #now also remove the individuals that don't have complete gps data (missing more than a minute) during the relevant time interval
  for(d in 1:nrow(labeled.intervals)){
    for(i in 1:n.inds){
      if(sum(is.na(allX[i,labeled.intervals$audio.start[d]:labeled.intervals$audio.end[d]])) > 60){
        inds.present[d,i] <- F
      }
    }
  }
  
  #if there are any dates where fewer than min.inds.present (default = 5) individuals were present, don't include them
  n.inds.present <- rowSums(inds.present)
  not.enough.inds <- which(n.inds.present < min.inds.present)
  if(length(not.enough.inds) > 0){
    print(paste('The interval for date', dates[not.enough.inds], 'contains fewer than the required number of individuals - removing'))
    dates <- dates[-not.enough.inds]
    inds.present <- inds.present[-not.enough.inds,]
    labeled.intervals <- labeled.intervals[-not.enough.inds,]
  }
  
  #filter to only focal calls - not needed with new pre-filtered data
  #calls.all <- calls.all[which(calls.all$focalType == 'F'),]
  
  #filter to only calls you want to analyze - currently close calls (cc) including cc hybrids or short notes (sn)
  calls.include <- calls.all[grep(callType,calls.all$callType),]
  
  #filter to only include calls within the recording intervals specified by date.intervals
  include.idxs <- c()
  for(date.idx in 1:length(dates)){
    include.idxs <- c(include.idxs, labeled.intervals$audio.start[date.idx]:labeled.intervals$audio.end[date.idx])
  }
  calls.include <- calls.include[which(calls.include$time.idx %in% include.idxs),]
  
  #-----------COMPUTE K STAT IN REAL VS RANDOMIZED DATA -----
  
  K.data <- matrix(NA, nrow = length(dist.windows), ncol = length(time.windows))
  print('computing K for real data')
  timestamp()
  out <- compute_Knox_indexes(calls.include, allX, allY, dist.windows, time.windows)
  K.data <- out$K
  tot.pairs.data <- out$tot.pairs
  
  #Randomizations - randomizing individuals on each day, only within individuals present on that day
  print('computing K for randomizations')
  timestamp()
  K.rand <- array(NA, dim = c(length(dist.windows)-1, length(time.windows)-1, n.rands))
  tot.pairs.rand <- rep(NA, n.rands)
  for(n in 1:n.rands){
    print(paste0('randomization ', n, '/', n.rands))
    allX_rand <- allX
    allY_rand <- allY
    
    for(date.idx in 1:length(dates)){
      inds.pres.date <- which(inds.present[date.idx,])
      n.matches <- 10
      while(n.matches > 0){
        permute.idxs <- sample(inds.pres.date, replace=F)
        n.matches <- sum(permute.idxs == inds.pres.date)
      }
      
      #replace GPS data with permuted data, for that date interval
      for(ind in 1:length(inds.pres.date)){
        time.idxs.date <- labeled.intervals$audio.start[date.idx]:labeled.intervals$audio.end[date.idx]
        allX_rand[inds.pres.date[ind], time.idxs.date] <- allX[permute.idxs[ind], time.idxs.date]
        allY_rand[inds.pres.date[ind], time.idxs.date] <- allY[permute.idxs[ind], time.idxs.date]
      }
    }
    
    out <- compute_Knox_indexes(calls.include, allX_rand, allY_rand, dist.windows, time.windows)
    K.rand[,,n] <- out$K
    tot.pairs.rand[n] <- out$tot.pairs
    timestamp()
  }
  
  #output file name
  outfile.name <- paste0(savedir,callType,'_clustering_',session,'_anon.RData')
  if(testflag){
    outfile.name <- paste0(savedir,callType,'_clustering_',session,'_test.RData')
  }
  
  #save
  save(list=c('session','dist.windows','time.windows','K.data','K.rand','tot.pairs.data','tot.pairs.rand','n.rands'), file = outfile.name)

}
    
    
