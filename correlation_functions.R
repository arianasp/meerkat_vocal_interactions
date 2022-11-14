#Correlation functions - are call rates correlated over space?

#Here we compute the correlation function (defined similar to Cavagna et al 2010)
#of call rate as a function of distance

#We first compute the difference from average call rate (at a given moment)
#across all individuals in the group, then measure the difference between
#an individual's call rate and this average call rate. This is what we use to
#compute the correlation function as a function of the radius r.

#--------------PARAMETERS----------

#call rate time window
call.rate.time.window <- 60

#distance bins 
dist.bins <- seq(0,30,3)

#list of sessions to use
sessions <- c('HM2017', 'HM2019', 'L2019')

#directories where data is stored
audiodir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/acoustic/'
gpsdir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/movement/'

#path to Baptiste's useful functions library
useful.fun.path <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/scripts/useful_functions.R'

#------------SETUP-------------
source(useful.fun.path)

#------------LOAD DATA AND COMPUTE C(r)-------
sess.idx <- 3

session <- sessions[sess.idx]

#load audio data
calls.all <- read.csv(paste0(audiodir,session, '_ALL_CALLS_SYNCHED.csv'), header=T, sep='\t', stringsAsFactors=F)

#loading spatial data and making names easier
longNames <- c(load(paste0(gpsdir,session,"_COORDINATES_all_sessions_level1.RData")))
shortNames <- simplifyNames(pat=paste(session,"_",sep=""))

#get number of individuals
n.inds <- nrow(indInfo)

#get number of time steps
n.times <- ncol(allX)

#call rates - close calls
windowCloseCalls <- t(apply(callTimeLine,1,FUN=function(f)(slidingWindow(f,function(g){
  if(all(is.na(g))){
    return(NA)
  }else{
    length(which(grepl("cc",g,ignore.case = T)))
  }
},window=call.rate.time.window,step=1,windowTime = "end"))))
call.rates <- windowCloseCalls / call.rate.time.window * 60 #convert to calls per minute

#compute dyadic distances
dyad.dists <- array(NA, dim = c(n.inds, n.inds, n.times))
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    dyad.dists[i,j,] <- sqrt((allX[i,] - allX[j,])^2 + (allY[i,] - allY[j,])^2)
  }
}
mean.dyad.dists <- apply(dyad.dists, 3, mean, na.rm=T)

#compute deviation from group average call rate
cc.avg <- apply(call.rates, 2, mean, na.rm=T)
cc.avg.mat <- matrix(rep(cc.avg, each = n.inds), nrow = n.inds, ncol = n.times)
cc.devs <- call.rates - cc.avg.mat

#compute the product of the deviations from average call rate for each pair
cc.prods <- array(NA, dim = c(n.inds, n.inds, n.times))
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    cc.prods[i,j,] <- cc.devs[i,] * cc.devs[j,]
  }
}

#for each bin, compute the average cc.prod
cc.cor <- counts <- matrix(NA, nrow = length(dist.bins)-1, ncol = n.times)
for(t in 1:n.times){
  dyad.dists.t <- dyad.dists[,,t]
  cc.prods.t <- cc.prods[,,t]
  for(d in 1:(length(dist.bins)-1)){
    idxs <- which(dyad.dists.t >= dist.bins[d] & dyad.dists.t < dist.bins[d+1])
    counts[d,t] <- length(idxs)
    cc.cor[d,t] <- mean(cc.prods.t[idxs], na.rm=T)
  }
}

#--------Mean dyadic distance vs average cc rate-----
mean.dyad.dist.bins <- seq(quantile(mean.dyad.dists,0.025, na.rm=T), quantile(mean.dyad.dists,0.975, na.rm=T), length.out = 10)
cc.rate.by.mean.dyad.dist <- rep(NA, length(mean.dyad.dist.bins)-1)
for(i in 1:(length(mean.dyad.dist.bins)-1)){
  idxs <- which(mean.dyad.dists >= mean.dyad.dist.bins[i] & mean.dyad.dists < mean.dyad.dist.bins[i+1])
  cc.rate.by.mean.dyad.dist[i] <- mean(cc.avg[idxs], na.rm=T)
}
plot(mean.dyad.dist.bins[2:length(mean.dyad.dist.bins)], cc.rate.by.mean.dyad.dist, pch = 19)

# quartz()
# plot(NULL, xlim = range(dist.bins), ylim = c(-20,20))
# for(i in 1:n.times){
#   lines(dist.bins[2:length(dist.bins)],cc.cor[,i], lwd = 0.01)
# }

#This method isn't working very well, probably because it is averaging over
#many time steps with varying distances among individuals. Need another
#method to measure how "clumped" calls are relative to one another and
#if it's more than expected based on where meerkats are.

#Idea: modification of "Knox metric".

#Compute the total number of pairs of events that are within a distance dr
#and a time window dt of one another. Divide this by the total number of 
#pairs of events that are within the time window dt of one another to get
#a value between 0 and 1. 0 = no events within distance window. 1 = all events
#within distance window.

#Now construct a null model by permuting call data (assign call data from
#individual i to individual j, ensuring i != j)

#Lastly, can vary both the time window dt and the radius dr to see how far
#and for how long call clusters persist.

