#Spatiotemporal call clustering - do calls cluster together in space?
#This code runs the call spatiotemporal clustering analysis described in the paper 
#"Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups"



#How does it work?
#It uses a modification of the "Knox metric".

#Compute the total number of pairs of events that are within a distance dr
#and a time window dt of one another. Divide this by the total number of 
#pairs of events that are within the time window dt of one another to get
#a value between 0 and 1. 0 = no events within distance window. 1 = all events
#within distance window.

#Now construct a null model by permuting call data (assign call data from
#individual i to individual j, ensuring i != j)

#Lastly, can vary both the time window dt and the radius dr to see how far
#and for how long call clusters persist.

#TODO: There are more pairs of calls in the data than in the randomizations...
#this is because sometimes individuals get their calls randomized to untracked individuals
#Possible solutions: normalize by total number of possible pairs somehow
#OR change the randomization somehow to only pair w/ inds tracked at that time (maybe by day?) <-- went with this one
#OR throw out randomization entirely and just look at raw values

#----------LIBRARIES--------
library(gplots)
library(fields)
library(viridis)

#--------------PARAMETERS----------

#time windows (sec)
time.windows <- c(1,5,10,30,60,120,300,600,1200,1800,3600)

#distance windows (m)
dist.windows <- c(2,4,6,8,10,15,20,25,30,40,60,80,100)

#number of randomizations
n.rands <- 100

#call type
callType <- 'cc'

#list of sessions to use
sessions <- c('HM2017', 'HM2019', 'L2019')

#directories where data is stored
audiodir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/acoustic/'
gpsdir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/movement/'
outdir <- '~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/data/call_clustering/'
plotdir <- paste0('~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/plots/',callType,'_call_clustering/')

#path to Baptiste's useful functions library
useful.fun.path <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/scripts/useful_functions.R'

#load pre-computed data?
load.precomputed.data <- T

#generate plots?
generate.plots <- T

#-------- FUNCTIONS ------

#Compute the modified Knox index (see description above)
#INPUTS:
# calls.include: call data frame
# allX, allY: matrices [n.inds x n.times] of x and y positions respectively
# dist.window: distance window used by the metric (m)
# time.window: time window used by the metric (sec)
#OUTPUTS:
# out: a list containing
#   out$dist.window: distance window (same as input)
#   out$time.window: time window (same as input)
#   out$num: the numerator of K, i.e. the total number of calls within a distance window dist.window and a time window time.window
#   out$denom: the denominator of K, i.e. the total number of calls within a time window time.window
#   out$K: the modified Knox K metric (num / denom)
#NOTE: In the calculations below, we actually don't use the raw K metric but rather take partial (numerical) derivatives using
#the num and denom outputs
compute_Knox_index <- function(calls.include, allX, allY, dist.window, time.window){
  
  #get x and y positions at the times of calls
  x_call <- allX[cbind(calls.include$ind.idx, calls.include$time.idx)]
  y_call <- allY[cbind(calls.include$ind.idx, calls.include$time.idx)]
  
  #construct a spatial distance matrix (actually formatted as a vector) for the dist between every pair of calls
  spatial.dist <- dist(cbind(x_call, y_call))
  
  #construct a temporal distance matrix for the temporal distance between every pair of calls
  temporal.dist <- dist(as.numeric(as.POSIXlt(calls.include$t0GPS_UTC)))
  
  #construct matrix to determine if the caller was the same caller (1) or a different (0)
  same.ind <- dist(calls.include$ind.idx) == 0
  
  #replace all instances of the same caller with NA
  spatial.dist[same.ind] <- NA
  temporal.dist[same.ind] <- NA
  
  #get the numerator of K,
  #this is the total number of calls within a distance window dist.window and a time window time.window
  num <- sum(spatial.dist <= dist.window & temporal.dist <= time.window, na.rm=T)
  
  #get the denominator (normalization factor) of K
  #this is the total number of calls within a time window time.window
  denom <- sum(temporal.dist <= time.window, na.rm=T)
  
  #compute K
  K <- num / denom
  
  #list for output
  out <- list()
  out$dist.window <- dist.window
  out$time.window <- time.window
  out$num <- num
  out$denom <- denom
  out$K <- K
  
  return(out)
}

#------------SETUP-------------
source(useful.fun.path)


if(!load.precomputed.data){
  print('Running call clustering analysis')
  timestamp()
  for(sess.idx in 1:length(sessions)){
  
    #------------LOAD DATA -------
    session <- sessions[sess.idx]
    
    print(session)
    
    #loading spatial data and making names easier
    longNames <- c(load(paste0(gpsdir,session,"_COORDINATES_all_sessions_level1.RData")))
    shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
    
    #load audio data
    calls.all <- read.csv(paste0(audiodir,session, '_ALL_CALLS_SYNCHED.csv'), header=T, sep='\t', stringsAsFactors=F)
    
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
    dates <- unique(calls.all$date)
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
          start.times.ind <- min(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T)
        }
        if(length(end.idxs)>0){
          end.times.ind <- calls.date$time.idx[end.idxs]
        } else{
          #if there were no 'start' labels found, use the earliest time of any entry for that individual
          end.times.ind <- max(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T)
        }
        
        #if there is more than one start label, use the latest one
        start.time.ind <- max(start.times.ind)
        
        #if there is a more than one end label, use the earliest one
        end.time.ind <- min(end.times.ind)
        
        #if the start or end time index is NA, this means it extends outside the bounds of the GPS data - replace it with the first / last non-NA time index point
        start.time.ind <- min(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T)
        end.time.ind <- max(calls.date$time.idx[which(calls.date$ind.idx==ind)], na.rm=T)
        
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
    empty.intervals <- which(labeled.intervals$audio.start >= labeled.intervals$audio.end)
    if(length(empty.intervals)>0){
      print(paste('The interval for date', dates[empty.intervals], 'ends before it starts - removing'))
      labeled.intervals <- labeled.intervals[-empty.intervals,]
      inds.present <- inds.present[-empty.intervals,]
      dates <- dates[-empty.intervals]
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
    for(r in 1:length(dist.windows)){
      for(t in 1:length(time.windows)){
        
        #compute K for the real data
        out <- compute_Knox_index(calls.include, allX, allY, dist.windows[r], time.windows[t])
        K.data[r,t] <- out$K
        num.data[r,t] <- out$num
        denom.data[r,t] <- out$denom
        
      }
    }
    
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
      
      #compute k index for each permutation
      for(r in 1:length(dist.windows)){
        for(t in 1:length(time.windows)){
          out <- compute_Knox_index(calls.include, allX_rand, allY_rand, dist.windows[r], time.windows[t])
          K.rand[r,t,n] <- out$K
          num.rand[r,t,n] <- out$num
          denom.rand[r,t,n] <- out$denom
        }
      }  
      
      
    }
    
    timestamp()
    
    #---------Calculate dK-----
    #dK gives the additional fraction within a certain spatial 'ring'
    #It is calculated by taking the total number of pairs of events within the time window and distance R
    #and subtracting the number of pairs of events within the time window and distance R_previous (the previous bin)
    #This value ('dnum' for 'difference in numberator') is then normalized by dividing by the total number of events in that time window (the denominator)
    # dnum.data <- num.data
    # dnum.data[2:dim(K.data)[1],] <- num.data[2:(dim(K.data)[1]),] - num.data[1:(dim(K.data)[1]-1),]
    # dK.data <- dnum.data / denom.data
    # 
    # dnum.rand <- num.rand
    # dnum.rand[2:dim(K.data)[1],,] <- num.rand[2:(dim(K.data)[1]),,] - num.rand[1:(dim(K.data)[1]-1),,]
    # dK.rand <- dnum.rand / denom.rand
    
    #--------Calculate dK / dA and dK / dt -----
    #dK / dA (r,t) is the spatial derivative of K, evaluated at time lag t
    #and distance r.
    #It gives the increase in number of pairs of calls captured when you 
    #epxand the spatial range from r to r + dr, divided by the area of the "ring"
    #between r and r + dr and the number of seconds elapsed t.
    #units of dK / dA are pairs of calls / sec / m^2. It is calculated by taking
    #the difference in rows of num.data and dividing this by the area (A) of 
    #the corresponding ring pi((r+dr)^2 - r^2) and by the time bin t.
    
    #dK / dt (r,t) is the temporal derivative of K, evaluated at time lag t
    #and radius r.
    #It gives the increase in pairs of calls for a given radius r when you
    #increase the time window from t to t + dt.
    #Units of dK / dt are pairs of calls / sec / m^2. It is calculated by 
    #taking the difference in columns of num.data and dividing this by dt 
    #and by the area of the current circle pi*r^2
    
    #rows = distances
    #cols = times
    #third dim = permutation number
    
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
    
    #Calculate the derivatives
    #dKdA i.e. partial derivative with respect to area (distance)
    dKdA.data <- (rbind(num.data[1,],num.data[2:nrow(num.data),] - num.data[1:(nrow(num.data)-1),]))/
      (dA.mat * time.windows.mat)
    dKdA.rand <- array(NA, dim = dim(num.rand))
    dKdA.rand[1,,] <- num.rand[1,,] / 
      (dA.array[1,,] * time.windows.array[1,,])
    dKdA.rand[2:nrow(num.data),,] <- (num.rand[2:nrow(num.data),,] - num.rand[1:(nrow(num.data)-1),,])/
      (dA.array[2:nrow(num.data),,] * time.windows.array[2:nrow(num.data),,])
    
    #dKdt i.e. partial derivative with respect to time
    dKdt.data <- (cbind(num.data[,1],num.data[,2:ncol(num.data)] - num.data[,1:(ncol(num.data)-1)]))/
      (circle.areas.mat * dt.mat)
    dKdt.rand <- array(NA, dim = dim(num.rand))
    dKdt.rand[,1,] <- num.rand[,1,] / 
      (circle.areas.array[,1,] * time.windows.array[,1,])
    dKdt.rand[,2:ncol(num.data),] <- (num.rand[,2:ncol(num.data),] - num.rand[,1:(ncol(num.data)-1),])/
      (dA.array[,2:ncol(num.data),] * time.windows.array[,2:ncol(num.data),])
    
    save(list=c('session','dist.windows','time.windows','dKdA.data','dKdA.rand','dKdt.data','dKdt.rand','num.data','num.rand','denom.data','denom.rand','n.rands'), file = paste0(outdir,callType,'_clustering_',session,'.RData'))
  
  }
}

#--------------PLOTTING---------

if(generate.plots){

  for(sess.idx in 1:length(sessions)){
  
    session <- sessions[sess.idx]
    
    setwd(outdir)
    load(paste0(outdir,callType, '_clustering_',session,'.RData'))
    
    setwd(plotdir)
    
    #NEW MATRICES - later move this up to the previous section, but here for now for plotting purposes
    
    #---This stuff will not be needed because it's computed above anyway
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
    
    #----end of stuff that won't be needed
    
    #---new stuff that will be needed to move up
    #calculate dt and store in a matrix / array
    dt.mat <- cbind(time.windows.mat[,1], time.windows.mat[,2:ncol(time.windows.mat)] - time.windows.mat[,1:(ncol(time.windows.mat)-1)])
    dt.array <- array(rep(dt.mat, n.rands), dim = c(length(dist.windows), length(time.windows), n.rands))
    #get new matrix / array, dKdAdt.data and dKdAdt.rand
    #this is the spatial and temporal derivative, in other words, 
    #each cell represents the density of pairs of calls (per area and per second)
    #that were given within a ring at R to R + dr and within the time window t to t + dt
    #It differs from dkdA in that the time windows only cover t to t + dt as opposed to covering all of t + dt
    #so it is more of any instantaneous temporal estimate of the clustering at that spatial scale as opposed to
    #"smeared" over a whole large time window
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
    
    #then normalize
    dKdAdt.data <- coldiffs.data / (dA.mat * dt.mat)
    dKdAdt.rand <- coldiffs.rand / (dA.array * dt.array)
    
    #---end new stuff
    
    #dKdAdt vs distance (lines = time)
    png(filename = paste0('dKdAdt_vs_dist_',session,'.png'), units = 'px', height = 800, width = 800)
    par(mfrow=c(1,1), mar = c(6,6,2,1))
    cols <- rev(viridis(length(time.windows)))
    plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdAdt.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (pairs / m^2 / sec)', log = 'x')
    for(t in 1:length(time.windows)){
      lines(dist.windows, dKdAdt.data[,t], lwd=3, col = cols[t], pch = 19)
    }
    legend('topright',legend = paste(time.windows,'sec'), col = cols, lty = 1, lwd = 3)
    dev.off()
    
    #dKdA vs distance (panels = time)
    png(filename = paste0('dKdA_vs_dist_comparenull',session,'.png'), units = 'px', height = 600, width = 2200)
    par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdA.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK / area)', main = paste('dt =', time.windows[t],'sec'), log = 'x')
      abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
      for(n in 1:n.rands){
        lines(dist.windows, dKdA.rand[,t,n], lwd = 0.5, col = 'black')
        #points(dist.windows, dK.rand[,t,n], cex = 0.5, col = 'black')
      }
      lines(dist.windows, dKdA.data[,t], lwd=2, col = 'red', pch = 19)
      #points(dist.windows, dK.data[,t], cex = 2, col = 'red', pch = 19)
    }
    dev.off()
    
    #dKdt vs distance (panels = time)
    png(filename = paste0('dKdt_vs_dist_comparenull',session,'.png'), units = 'px', height = 600, width = 2200)
    par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdt.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK / area)', main = paste('dt =', time.windows[t],'sec'), log = 'x')
      abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
      for(n in 1:n.rands){
        lines(dist.windows, dKdt.rand[,t,n], lwd = 0.5, col = 'black')
        #points(dist.windows, dK.rand[,t,n], cex = 0.5, col = 'black')
      }
      lines(dist.windows, dKdt.data[,t], lwd=2, col = 'red', pch = 19)
      #points(dist.windows, dK.data[,t], cex = 2, col = 'red', pch = 19)
    }
    dev.off()
    
    #dK/(dAdt) vs distance (panels = time)
    png(filename = paste0('dKdAdt_vs_dist_comparenull',session,'.png'), units = 'px', height = 600, width = 2200)
    par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdAdt.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK / area)', main = paste('dt =', time.windows[t],'sec'), log = 'x')
      abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
      for(n in 1:n.rands){
        lines(dist.windows, dKdAdt.rand[,t,n], lwd = 0.5, col = 'black')
        #points(dist.windows, dK.rand[,t,n], cex = 0.5, col = 'black')
      }
      lines(dist.windows, dKdAdt.data[,t], lwd=2, col = 'red', pch = 19)
      #points(dist.windows, dK.data[,t], cex = 2, col = 'red', pch = 19)
    }
    dev.off()
    
    #dK/(dt) vs distance (panels = time)
    png(filename = paste0('dKdt_vs_dist_comparenull',session,'.png'), units = 'px', height = 600, width = 2200)
    par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdA.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK / area)', main = paste('dt =', time.windows[t],'sec'), log = 'x')
      abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
      for(n in 1:n.rands){
        lines(dist.windows, dKdt.rand[,t,n], lwd = 0.5, col = 'black')
        #points(dist.windows, dK.rand[,t,n], cex = 0.5, col = 'black')
      }
      lines(dist.windows, dKdt.data[,t], lwd=2, col = 'red', pch = 19)
      #points(dist.windows, dK.data[,t], cex = 2, col = 'red', pch = 19)
    }
    dev.off()
    
    #dKdAdt vs time and distance, log(data / mean(null))
    png(filename = paste0('dKdAdt_vs_dist_and_time_comparenull_logratio',session,'.png'), units = 'px', height = 800, width = 800)
    par(mfrow=c(1,1), cex.main = 2, cex.lab=2, mar = c(5,5,1,1), cex.axis = 1.5)
    dKdAdt.rand.mean <- apply(dKdAdt.rand, c(1,2), mean)
    ratio <- log(dKdAdt.data / dKdAdt.rand.mean, base = 10)
    image.plot(ratio, zlim = c(-max(ratio,na.rm=T),max(ratio,na.rm=T)), col = bluered(256), xlab = 'Distance (m)', ylab = 'Time (sec)', xaxt = 'n', yaxt = 'n')
    axis(side = 1, at = seq(0,1,length.out = length(dist.windows)), labels = dist.windows)
    axis(side = 2, at = seq(0,1,length.out = length(time.windows)), labels = time.windows)
    dev.off()
  }
}
    
    
    
