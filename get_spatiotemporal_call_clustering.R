#Spatiotemporal call clustering - do calls cluster together in space?

#This code runs the call spatiotemporal clustering analysis described in the paper 
#"Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups"
#It starts from the raw GPS and call label data and outputs processed data that are used in a later plotting 
#script (plot_spatiotemporal_call_clustering.R)

#How does it work?
#First, we define a clustering metric.
#Here we use modification of the "Knox statistic" (actually kind of numerical derivative of the Knox statistic, described below)
#Basically, as a metric of clustering, we compute the number of pairs of calls that 
#fall within a distance range r to r + dr and a temporal range t to t + dt. We then noramlize this value by the
#spatial and temporal bin size to get a kind of "density" estimate. We divide by dt, and also by the area of the
#"ring" between the circle with radius r and the circle with radius r + dr. So for the spatial part, we divide
#by the value pi * ((r+dr)^2 - (r)^2). Why a ring rather than just dividing by dr? Because there are more ways
#for pairs of points to be larger distances apart, so this needs to be normalized by an area instead of by a distance.
#In this case the relevant area is the ring described above. 
#So, overall, the metric is defined as
# K(r,t) = number of pairs of calls / [pi * ((r+dr)^2 - (r)^2)] / [dt]
#Its units are pairs of calls per m^2 per second

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
audiodir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 
gpsdir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 

#directory where you would like to store the processed data (can be the same directory too if you like)
savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/'

#call type
callType <- 'cc'

#test on a smaller number of bins and with only one null model?
testflag <- F

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------

#time windows (sec)
time.windows <- c(1,3,10,30,90,180,600,1800,5400)

#distance windows (m)
dist.windows <- c(1,2,5,10,15,20,30,50,100)

#number of randomizations
n.rands <- 100

#list of sessions to use
sessions <- c('HM2019','HM2017', 'L2019')

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
#OUTPUTS:
# out: a list containing
#   out$dist.window: distance window (same as input)
#   out$time.window: time window (same as input)
#   out$num: matrix of the numerators of K for each distance (row) / time window (column), i.e. the total number of calls within a distance window dist.window and a time window time.window
#   out$denom: matrix of the denominators of K, i.e. the total number of calls within a time window time.window
#   out$K: matrix of the modified Knox K metric (num / denom)
#NOTE: In the calculations below, we actually don't use the raw K metric but rather take partial (numerical) derivatives using
#the num and denom outputs
compute_Knox_indexes <- function(calls.include, allX, allY, dist.windows, time.windows){
  
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
  
  #construct matrix to determine if the caller was the same caller (1) or a different (0)
  diff.inds <- dist(calls.include$ind.idx) != 0
  
  #remove same caller data
  spatial.dist <- spatial.dist[diff.inds]
  temporal.dist <- temporal.dist[diff.inds]
  
  #matrices to store output
  Ks <- nums <- denoms <- matrix(NA, nrow = length(dist.windows), ncol = length(time.windows))
  for(r in 1:length(dist.windows)){
    for(t in 1:length(time.windows)){

      #get the numerator of K,
      #this is the total number of calls within a distance window dist.window and a time window time.window
      nums[r,t] <- sum(spatial.dist <= dist.windows[r] & temporal.dist <= time.windows[t], na.rm=T)
      
      #get the denominator (normalization factor) of K
      #this is the total number of calls within a time window time.window
      denoms[r,t] <- sum(temporal.dist <= time.windows[t], na.rm=T)
      
      #compute K
      Ks[r,t] <- nums[r,t] / denoms[r,t]
      
    }
  }
  
  #list for output
  out <- list()
  out$dist.windows <- dist.windows
  out$time.windows <- time.windows
  out$nums <- nums
  out$denoms <- denoms
  out$Ks <- Ks
  
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
calls.allsessions <- read.csv(paste0(audiodir,'all_calls_sync_resolved_with_oor_2022-12-04.csv' ), header=T, stringsAsFactors=F)

timestamp()
for(sess.idx in 1:length(sessions)){

  #------------LOAD DATA -------
  session <- sessions[sess.idx]
  
  print(session)
  
  #loading spatial data and making names easier
  longNames <- c(load(paste0(gpsdir,session,"_COORDINATES_all_sessions.RData")))
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
  #Note: This will not account for instances where there are gaps in labeling (but this should be rare, unlikely to systematically affect this analysis as it will be same in real vs randomized data)
  labeled.intervals <- data.frame(date = dates, 
                               audio.start = rep(NA,length(dates)),
                               audio.end = rep(NA, length(dates)))
  inds.present <- matrix(FALSE, nrow = nrow(labeled.intervals), ncol = n.inds)
  for(date.idx in 1:length(dates)){
    calls.date <- calls.all[which(calls.all$date == dates[date.idx]),]
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
  
  #if there are any dates where fewer than min.inds.present (default = 5) individuals were present, don't include them
  n.inds.present <- rowSums(inds.present)
  not.enough.inds <- which(n.inds.present < min.inds.present)
  if(length(not.enough.inds) > 0){
    print(paste('The interval for date', dates[not.enough.inds], 'contains fewer than the required number of individuals - removing'))
    dates <- dates[-not.enough.inds]
    inds.present <- inds.present[-not.enough.inds,]
    labeled.intervals <- labeled.intervals[-not.enough.inds,]
  }
  
  #filter to only focal calls
  calls <- calls.all[which(calls.all$focalType == 'F'),]
  
  #filter to only calls you want to analyze - currently close calls (cc) including cc hybrids or short notes (sn)
  calls.include <- calls[grep(callType,calls$callType),]
  
  #filter to only include calls within the recording intervals specified by date.intervals
  include.idxs <- c()
  for(date.idx in 1:length(dates)){
    include.idxs <- c(include.idxs, labeled.intervals$audio.start[date.idx]:labeled.intervals$audio.end[date.idx])
  }
  calls.include <- calls.include[which(calls.include$time.idx %in% include.idxs),]
  
  #-----------COMPUTE K STAT IN REAL VS RANDOMIZED DATA -----
  
  K.data <- num.data <- denom.data <- matrix(NA, nrow = length(dist.windows), ncol = length(time.windows))
  print('computing K for real data')
  timestamp()
  out <- compute_Knox_indexes(calls.include, allX, allY, dist.windows, time.windows)
  K.data <- out$K
  num.data <- out$num
  denom.data <- out$denom
  
  #Randomizations - randomizing individuals on each day, only within individuals present on that day
  print('computing K for randomizations')
  timestamp()
  K.rand <- num.rand <- denom.rand <- array(NA, dim = c(length(dist.windows), length(time.windows), n.rands))
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
    num.rand[,,n] <- out$num
    denom.rand[,,n] <- out$denom
    timestamp()
  }
  
  #calculate circle areas and convert to a matrix (for data) and array (for randomizations)
  circle.areas <- pi*dist.windows^2
  circle.areas.mat <- matrix(rep(circle.areas, length(time.windows)), nrow = length(dist.windows), ncol = length(time.windows))
  circle.areas.array <- array(rep(circle.areas.mat, n.rands), dim = c(length(dist.windows), length(time.windows), n.rands))
  
  #construct a time matrix
  time.windows.mat <- matrix(rep(time.windows, each = length(dist.windows)), nrow = length(dist.windows), ncol = length(time.windows))
  time.windows.array <- array(rep(time.windows.mat, n.rands), dim = c(length(dist.windows), length(time.windows), n.rands))
  
  #calculate differences between circle areas and store in a matrix / array
  dA.mat <- rbind(circle.areas.mat[1,], (circle.areas.mat[2:nrow(circle.areas.mat),] - circle.areas.mat[1:(nrow(circle.areas.mat)-1)]))
  dA.array <- array(rep(dA.mat, n.rands), dim = c(length(dist.windows), length(time.windows), n.rands))
  
  #calculate dt and store in a matrix / array
  dt.mat <- cbind(time.windows.mat[,1], time.windows.mat[,2:ncol(time.windows.mat)] - time.windows.mat[,1:(ncol(time.windows.mat)-1)])
  dt.array <- array(rep(dt.mat, n.rands), dim = c(length(dist.windows), length(time.windows), n.rands))
  
  #get new matrix / array, dKdAdt.data and dKdAdt.rand
  #this is the spatial and temporal derivative, in other words, 
  #each cell represents the density of pairs of calls (per area and per second)
  #that were given within a ring at R to R + dr and within the time window t to t + dt
  #rows = distances
  #cols = times
  #third dim = permutation number
  
  #take the row (distance) differences in number of call pairs
  rowdiffs.data <- (rbind(num.data[1,],num.data[2:nrow(num.data),] - num.data[1:(nrow(num.data)-1),]))
  rowdiffs.rand <- array(NA, dim = dim(num.rand))
  rowdiffs.rand[1,,] <- num.rand[1,,]
  rowdiffs.rand[2:nrow(num.data),,] <- (num.rand[2:nrow(num.data),,] - num.rand[1:(nrow(num.data)-1),,])
  
  #then take the column (time) differences in number of call pairs
  coldiffs.data <- (cbind(rowdiffs.data[,1],rowdiffs.data[,2:ncol(rowdiffs.data)] - rowdiffs.data[,1:(ncol(rowdiffs.data)-1)]))
  coldiffs.rand <- array(NA, dim = dim(rowdiffs.rand))
  coldiffs.rand[,1,] <- rowdiffs.rand[,1,]
  coldiffs.rand[,2:ncol(rowdiffs.rand),] <- (rowdiffs.rand[,2:ncol(rowdiffs.rand),] - rowdiffs.rand[,1:(ncol(rowdiffs.rand)-1),])
  
  #then normalize by area and time window
  dKdAdt.data <- coldiffs.data / (dA.mat * dt.mat)
  dKdAdt.rand <- coldiffs.rand / (dA.array * dt.array)
  
  #output file name
  outfile.name <- paste0(savedir,callType,'_clustering_',session,'.RData')
  if(testflag){
    outfile.name <- paste0(savedir,callType,'_clustering_',session,'_test.RData')
  }
  
  #save
  save(list=c('session','dist.windows','time.windows','num.data','num.rand','denom.data','denom.rand','n.rands','dKdAdt.data','dKdAdt.rand'), file = outfile.name)

}
    
    
