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
mean.dyad.dists <- apply(dyad.dists, 2, mean, na.rm=T)

#compute deviation from group average call rate
cc.avg <- apply(call.rates, 2, mean, na.rm=T)
cc.avg.mat <- matrix(rep(cc.avg, each = n.inds), nrow = n.inds, ncol = n.times)
cc.devs <- call.rates - cc.avg.mat

#compute the product of the call rates for each pair
cc.prods <- array(NA, dim = c(n.inds, n.inds, n.times))
for(i in 1:(n.inds-1)){
  for(j in (i+1):n.inds){
    cc.prods[i,j,] <- cc.devs[i,] * cc.devs[j,]
  }
}

#for each bin, compute the average cc.prod
cc.cor <- counts <- rep(NA, length(dist.bins)-1)
for(d in 1:(length(dist.bins)-1)){
  idxs <- which(dyad.dists >= dist.bins[d] & dyad.dists < dist.bins[d+1])
  counts[d] <- length(idxs)
  cc.cor[d] <- mean(cc.prods[idxs], na.rm=T)
}

plot(dist.bins[2:length(dist.bins)], cc.cor, xlab = 'r (m)', ylab = 'C(r)', pch = 19, cex = counts / max(counts)*2)
abline(h=0)
