#Spatiotemporal call clustering - do calls cluster together in space?

#Modification of "Knox metric".

#Compute the total number of pairs of events that are within a distance dr
#and a time window dt of one another. Divide this by the total number of 
#pairs of events that are within the time window dt of one another to get
#a value between 0 and 1. 0 = no events within distance window. 1 = all events
#within distance window.

#Now construct a null model by permuting call data (assign call data from
#individual i to individual j, ensuring i != j)

#Lastly, can vary both the time window dt and the radius dr to see how far
#and for how long call clusters persist.

#--------------PARAMETERS----------

#time windows (sec)
time.windows <- c(1,5,10,30,60,120,300,600,1200)

#distance windows (m)
dist.windows <- c(2,4,6,8,10,15,20,25,30,40,60,80,100)

#number of randomizations
n.rands <- 10

#list of sessions to use
sessions <- c('HM2017', 'HM2019', 'L2019')

#directories where data is stored
audiodir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/acoustic/'
gpsdir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/movement/'
outdir <- '~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/data/call_clustering/'

#path to Baptiste's useful functions library
useful.fun.path <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/scripts/useful_functions.R'

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
# K: the modified Knox index (see above)
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
  
  #filter to only focal calls
  calls <- calls.all[which(calls.all$focalType == 'F'),]
  
  #filter to only calls you want to analyze
  calls.include <- calls[grep('cc',calls$callType),]
  
  #get number of individuals
  n.inds <- nrow(indInfo)
  
  #get number of time steps
  n.times <- ncol(allX)
  
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
  
  #randomizations
  print('computing K for randomizations')
  timestamp()
  K.rand <- num.rand <- denom.rand <- array(NA, dim = c(length(dist.windows), length(time.windows), n.rands))
  for(n in 1:n.rands){
    print(paste('n =', n, '/', n.rands))
    #permute indexes
    n.matches <- 10
    while(n.matches > 0){
      permute.idxs <- sample(n.inds, replace = F)
      n.matches <- sum(permute.idxs == seq(1:n.inds))
    }
    
    #compute k index for each permutation
    for(r in 1:length(dist.windows)){
      for(t in 1:length(time.windows)){
        out <- compute_Knox_index(calls.include, allX[permute.idxs,], allY[permute.idxs,], dist.windows[r], time.windows[t])
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
  dnum.data <- num.data
  dnum.data[2:dim(K.data)[1],] <- num.data[2:(dim(K.data)[1]),] - num.data[1:(dim(K.data)[1]-1),]
  dK.data <- dnum.data / denom.data
  
  dnum.rand <- num.rand
  dnum.rand[2:dim(K.data)[1],,] <- num.rand[2:(dim(K.data)[1]),,] - num.rand[1:(dim(K.data)[1]-1),,]
  dK.rand <- dnum.rand / denom.rand
  
  save(list=c('session','dist.windows','time.windows','dK.data','dK.rand','K.data','K.rand','num.data','num.rand','denom.data','denom.rand','n.rands'), file = paste0(outdir,'cc_clustering_',session,'.RData'))

}

#--------------PLOTTING---------

if(generate.plots){

  for(sess.idx in 1:length(sessions)){
  
    session <- sessions[sess.idx]
    
    setwd(outdir)
    load(paste0(outdir,'cc_clustering_',session,'.RData'))
    
    #Make a plot of data K vs randomizations
    quartz(height = 8, width = 12)
    par(mfrow=c(1,length(time.windows)))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,1), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance window (m)', ylab = 'K')
      for(n in 1:n.rands){
        points(dist.windows, K.rand[,t,n], cex = 0.5, col = 'black')
      }
      points(dist.windows, K.data[,t], cex = 2, col = 'red')
    }
    
    #Make a plot of data dK vs randomizations
    quartz(height = 6, width = 22)
    par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
    for(t in 1:length(time.windows)){
      plot(NULL, xlim = range(dist.windows), ylim = c(0,max(c(c(dK.data),c(dK.rand)),na.rm=T)*1.1), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK)', main = paste('dt =', time.windows[t],'sec'))
      abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
      for(n in 1:n.rands){
        lines(dist.windows, dK.rand[,t,n], lwd = 0.5, col = 'black')
        #points(dist.windows, dK.rand[,t,n], cex = 0.5, col = 'black')
      }
      lines(dist.windows, dK.data[,t], lwd=2, col = 'red', pch = 19)
      #points(dist.windows, dK.data[,t], cex = 2, col = 'red', pch = 19)
    }
  }
}
    
    
    
