#Exploring the patterns of calling at the group level

#TODO:
#Update script to use the most recent data

#----------------PARAMETERS--------------------

#whether to save output
save.output <- T

#directories where data is stored
audiodir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/acoustic/'
gpsdir <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/data/movement/'

#directory where code is stored for this project
codedir <- '~/Dropbox/code_ari/meerkat_vocal_interactions'

#path to Baptiste's useful functions library
useful.fun.path <- '~/Dropbox/code_ari/Meerkat_leadership_hackathon/scripts/useful_functions.R'

#directory of where to save results
savedir <- '~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/data/call_response'

#filename of general meerkat functions
general.funcs.filename <- '~/Dropbox/code_ari/move_comm_analysis/general/meerkat_functions.R'

#window over which to compute call rates
call.rate.time.window <- 60 #seconds

#spatial discretization step for centroid (relative positioning computation)
spatial.discretization.step <- 10

#spatial bins
spatial.bins <- seq(-30,30,5)

#normalized spatial bins
norm.spatial.bins <- seq(-2,2,length.out = 11)

#list of sessions to use
sessions <- c('HM2017', 'HM2019', 'L2019')

#call rate bins
cc.rate.bins <- c(0,1,2,3,1000)
move.rate.bins <- c(0,1,1000)
sn.rate.bins <- c(0,1,1000)

#nearest neighbor distance bins
nn.dist.bins <- c(0,2,5,10,10000)

#----------------SETUP-----------------

library(fields)
library(colorRamps)

#set computer time zone to UTC, to avoid stupid time zone issues
Sys.setenv(TZ='UTC')

#source useful functions
source(useful.fun.path)
source(general.funcs.filename)
#source('~/Dropbox/code_ari/Meerkat_leadership_hackathon/scripts/leadership_plotting_functions.R')

#---------------LOAD DATA---------------

#specify spatial bins
bins <- norm.spatial.bins #bins to use
mids <- (bins[1:(length(bins)-1)] + bins[2:length(bins)]) / 2

#initialize vectors to hold output

#call rates vs. spatial position
cc.x.means <- cc.y.means <- matrix(NA, nrow = length(sessions), ncol = length(mids))
mov.x.means <- mov.y.means <- matrix(NA, nrow = length(sessions), ncol = length(mids))
sn.x.means <- sn.y.means <- matrix(NA, nrow = length(sessions), ncol = length(mids))
x.counts <- y.counts <- matrix(NA, nrow = length(sessions), ncol = length(mids))

#correlations between call rates
cor.mats <- list()

#call rate vs nearest neighbor distance and nearest neighbor call rate
call.rate.vs.nn.mats <- list()

for(sess.idx in 1:length(sessions)){
  session <- sessions[sess.idx]

  timestamp()
  print('session:')
  print(session)

  #load audio data
  calls.all <- read.csv(paste0(audiodir,session, '_ALL_CALLS_SYNCHED.csv'), header=T, sep='\t', stringsAsFactors=F)
  
  #loading spatial data and making names easier
  longNames <- c(load(paste0(gpsdir,session,"_COORDINATES_all_sessions_level1.RData")))
  shortNames <- simplifyNames(pat=paste(session,"_",sep=""))
  
  #get number of individuals
  n.inds <- nrow(indInfo)
  
  #get number of time steps
  n.times <- ncol(allX)
  
  #--------------PROCESS-------------
  #call rates - move calls
  windowMoveCalls <- t(apply(callTimeLine,1,FUN=function(f)(slidingWindow(f,function(g){
    if(all(is.na(g))){
      return(NA)
    }else{
      length(which(grepl("mo",g,ignore.case = T) | grepl("ld",g,ignore.case = T)))
    }
  },window=call.rate.time.window,step=1,windowTime = "end"))))
  
  #call rates - close calls
  windowCloseCalls <- t(apply(callTimeLine,1,FUN=function(f)(slidingWindow(f,function(g){
    if(all(is.na(g))){
      return(NA)
    }else{
      length(which(grepl("cc",g,ignore.case = T)))
    }
  },window=call.rate.time.window,step=1,windowTime = "end"))))
  
  #call rates - short note calls
  windowShortNoteCalls <- t(apply(callTimeLine,1,FUN=function(f)(slidingWindow(f,function(g){
    if(all(is.na(g))){
      return(NA)
    }else{
      length(which(grepl("sn",g,ignore.case = T)))
    }
  },window=call.rate.time.window,step=1,windowTime = "end"))))
  
  #relative positioning
  spatialPast <- relativePos(allX,allY,step=spatial.discretization.step, discrSpatial=T,futur=F,timeline=timeLine,discrByCentroid=T,centroidSpeed = T, removeInd = F) #metrics based on past
  relX <- spatialPast$relative_ind_X
  relY <- spatialPast$relative_ind_Y
  sdX <- matrix(rep(apply(relX, 2, sd, na.rm=T), each = n.inds), nrow = n.inds, ncol = n.times)
  sdY <- matrix(rep(apply(relY, 2, sd, na.rm=T), each = n.inds), nrow = n.inds, ncol = n.times)
  normRelX <- relX / sdX
  normRelY <- relY / sdY
  
  #-------------ANALYSES-------------
  
  #----Call rate as a function of front-back (relative) distance
  X <- normRelX
  for(i in 1:length(mids)){
    curr <- which(X >= bins[i] & X < bins[i+1])
    cc.x.means[sess.idx,i] <- mean(windowCloseCalls[curr], na.rm=T)
    mov.x.means[sess.idx,i] <- mean(windowMoveCalls[curr], na.rm=T)
    sn.x.means[sess.idx,i] <- mean(windowShortNoteCalls[curr], na.rm=T)
    x.counts[sess.idx,i] <- sum(!is.na(windowCloseCalls[curr]))
  }
  
  #----Call rate as a function of left-right (relative) distance
  Y <- normRelY
  for(i in 1:length(mids)){
    curr <- which(Y >= bins[i] & Y < bins[i+1])
    cc.y.means[sess.idx,i] <- mean(windowCloseCalls[curr], na.rm=T)
    mov.y.means[sess.idx,i] <- mean(windowMoveCalls[curr], na.rm=T)
    sn.y.means[sess.idx,i] <- mean(windowShortNoteCalls[curr], na.rm=T)
    y.counts[sess.idx,i] <- sum(!is.na(windowCloseCalls[curr]))
  }
  
  #----Call rate correlations
  cor.cc.cc <- cor.cc.mov <- cor.cc.sn <- cor.mov.mov <- cor.mov.sn <- cor.sn.sn <- cor.counts <- matrix(NA, nrow = n.inds, ncol = n.inds)
  
  #cc vs cc
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowCloseCalls
  callMat2 <- windowCloseCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.cc.cc <- corMat
  
  #cc vs mov
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowCloseCalls
  callMat2 <- windowMoveCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.cc.mov <- corMat
  
  #cc vs sn
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowCloseCalls
  callMat2 <- windowShortNoteCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.cc.sn <- corMat

  #mov vs mov
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowMoveCalls
  callMat2 <- windowMoveCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.mov.mov <- corMat
    
  #mov vs sn
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowCloseCalls
  callMat2 <- windowCloseCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.mov.sn <- corMat
  
  #sn vs sn
  corMat <- countMat <- matrix(NA, nrow = n.inds, ncol = n.inds)
  callMat <- windowShortNoteCalls
  callMat2 <- windowShortNoteCalls
  for(i in 1:n.inds){
    for(j in 1:n.inds){
      npairs <- sum(!is.na(callMat[i,]) & !is.na(callMat2[j,]))
      countMat[i,j] <- npairs
      if(npairs > 100){
        corMat[i,j] <- cor(callMat[i,], callMat2[j,], use = 'complete', method = 'spearman')
      }
    }
  }
  cor.sn.sn <- corMat
  
  #counts
  cor.counts <- countMat #same across all call types
  
  cor.mats[[sess.idx]] <- list(cor.cc.cc = cor.cc.cc, cor.cc.mov = cor.cc.mov, cor.cc.sn = cor.cc.sn, 
                               cor.mov.mov = cor.mov.mov, cor.mov.sn = cor.mov.sn, 
                               cor.sn.sn = cor.sn.sn, 
                               cor.counts = cor.counts)
  
  #---Call rate as a function of nn distance and nn call rate
  
  #get dyadic distances
  dyad.dists <- array(NA, dim = c(nrow(allX),nrow(allX),ncol(allX)))
  for(ind1 in 1:(nrow(allX)-1)){
    for(ind2 in (ind1+1):nrow(allX)){
      dds <- sqrt((allX[ind1,] - allX[ind2,])^2 + (allY[ind1,] - allY[ind2,])^2)
      dyad.dists[ind1,ind2,] <- dds
      dyad.dists[ind2,ind1,] <- dds
    }
  }
  
  #get nearest neighbor distance and id
  nn.dists <- apply(dyad.dists,MARGIN = c(2,3), FUN = min, na.rm=T)
  nn.dists[which(is.infinite(nn.dists))] <- NA
  nn.ids <- apply(dyad.dists, MARGIN = c(2,3), FUN = function(x){return(which.min(x)[1])})
  
  #get nearest neighbor call rates
  nn.cc.rate <- nn.move.rate <- nn.sn.rate <- matrix(NA, nrow = nrow(allX), ncol = ncol(allX))
  for(ind in 1:nrow(allX)){
    nn.cc.rate[ind,] <- windowCloseCalls[cbind(nn.ids[ind,],1:ncol(allX))]
    nn.move.rate[ind,] <- windowMoveCalls[cbind(nn.ids[ind,],1:ncol(allX))]
    nn.sn.rate[ind,] <- windowShortNoteCalls[cbind(nn.ids[ind,],1:ncol(allX))]
  }
  
  #get mean close call rate as a function of nearest neighbor distance and call rate
  cc.rate.vs.nn.mat <- matrix(NA, nrow = length(cc.rate.bins)-1, ncol = length(nn.dist.bins)-1)
  for(rate in 1:(length(cc.rate.bins)-1)){
    for(dist in 1:(length(nn.dist.bins)-1)){
      idxs <- which(nn.dists >= nn.dist.bins[dist] & nn.dists < nn.dist.bins[dist+1] & nn.cc.rate >= cc.rate.bins[rate] & nn.cc.rate < cc.rate.bins[rate+1])
      if(length(idxs)>0){
        cc.rate.vs.nn.mat[rate, dist] <- mean(windowCloseCalls[idxs], na.rm=T)
      }
    }
  }
  
  #get mean move call rate as a function of nearest neighbor distance and call rate
  move.rate.vs.nn.mat <- matrix(NA, nrow = length(move.rate.bins)-1, ncol = length(nn.dist.bins)-1)
  for(rate in 1:(length(cc.rate.bins)-1)){
    for(dist in 1:(length(nn.dist.bins)-1)){
      idxs <- which(nn.dists >= nn.dist.bins[dist] & nn.dists < nn.dist.bins[dist+1] & nn.move.rate >= move.rate.bins[rate] & nn.move.rate < move.rate.bins[rate+1])
      if(length(idxs)>0){
        move.rate.vs.nn.mat[rate, dist] <- mean(windowMoveCalls[idxs], na.rm=T)
      }
    }
  }
  
  #get mean sn call rate as a function of nearest neighbor distance and call rate
  sn.rate.vs.nn.mat <- matrix(NA, nrow = length(sn.rate.bins)-1, ncol = length(nn.dist.bins)-1)
  for(rate in 1:(length(sn.rate.bins)-1)){
    for(dist in 1:(length(nn.dist.bins)-1)){
      idxs <- which(nn.dists >= nn.dist.bins[dist] & nn.dists < nn.dist.bins[dist+1] & nn.sn.rate >= sn.rate.bins[rate] & nn.sn.rate < sn.rate.bins[rate+1])
      if(length(idxs)>0){
        sn.rate.vs.nn.mat[rate, dist] <- mean(windowShortNoteCalls[idxs], na.rm=T)
      }
    }
  }
  
  call.rate.vs.nn.mats.sess <- list()
  call.rate.vs.nn.mats.sess$cc.rate.vs.nn.mat <- cc.rate.vs.nn.mat
  call.rate.vs.nn.mats.sess$move.rate.vs.nn.mat <- move.rate.vs.nn.mat
  call.rate.vs.nn.mats.sess$sn.rate.vs.nn.mat <- sn.rate.vs.nn.mat
  call.rate.vs.nn.mats.sess$cc.rate.bins <- cc.rate.bins
  call.rate.vs.nn.mats.sess$move.rate.bins <- move.rate.bins
  call.rate.vs.nn.mats.sess$sn.rate.bins <- sn.rate.bins
  call.rate.vs.nn.mats.sess$nn.dist.bins <- nn.dist.bins
   
  call.rate.vs.nn.mats[[sess.idx]] <- call.rate.vs.nn.mats.sess
  
}

call.position.data <- list(cc.x.means = cc.x.means, cc.y.means = cc.y.means, 
                           mov.x.means = mov.x.means, mov.y.means = mov.y.means, 
                           sn.x.means = sn.x.means, sn.y.means = sn.y.means,
                           x.counts = x.counts, y.counts = y.counts)

#----------PLOTS-------------
cols <- c('red','blue','darkgreen')

quartz(height = 12, width = 8)
par(mar=c(6,6,1,1), mfrow = c(3,2))

scalefactor <- max(c(x.counts,y.counts)) / 4

#cc front-back
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(cc.x.means,cc.y.means))), xlab = 'Front-back position (normalized)', ylab = 'Mean close call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, cc.x.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, cc.x.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = x.counts / scalefactor)
}

#cc left-right
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(cc.x.means,cc.y.means))), xlab = 'Left-right position (normalized)', ylab = 'Mean close call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, cc.y.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, cc.y.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = y.counts / scalefactor)
}

#mov front-back
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(mov.x.means,mov.y.means))), xlab = 'Front-back position (normalized)', ylab = 'Mean move call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, mov.x.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, mov.x.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = x.counts / scalefactor)
}
legend('topleft',legend = sessions, fill = cols, cex = 2)

#mov left-right
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(mov.x.means,mov.y.means))), xlab = 'Left-right position (normalized)', ylab = 'Mean move call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, mov.y.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, mov.y.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = y.counts / scalefactor)
}

#sn front-back
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(sn.x.means,sn.y.means))), xlab = 'Front-back position (normalized)', ylab = 'Mean short note call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, sn.x.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, sn.x.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = x.counts / scalefactor)
}

#sn left-right
plot(NULL, xlim = range(bins), ylim = c(0,max(cbind(sn.x.means,sn.y.means))), xlab = 'Left-right position (normalized)', ylab = 'Mean short note call rate (calls / min)',cex.lab=2, cex.axis = 1.5)
for(sess.idx in 1:length(sessions)){
  lines(mids, sn.y.means[sess.idx, ], col = cols[sess.idx],lwd=3)
  points(mids, sn.y.means[sess.idx, ], col = cols[sess.idx], pch = 19, cex = y.counts / scalefactor)
}


# #Are individuals correlated in their call rates?

pal <- colorRampPalette(c('red','white','blue'))
cols <- pal(256)

#CC vs all
quartz(width = 10, height = 10)
par(mfcol=c(length(sessions),3), mar = c(4,4,2,4))
for(i in 1:length(sessions)){
  non.na <- which(rowSums(!is.na(cor.mats[[i]]$cor.cc.cc)) >0)
  image.plot(cor.mats[[i]]$cor.cc.cc[non.na,non.na], col = cols, zlim = c(-1,1), main = 'CC vs CC')
  image.plot(cor.mats[[i]]$cor.cc.mov[non.na,non.na], col = cols, zlim = c(-1,1), main = 'CC vs MOV')
  image.plot(cor.mats[[i]]$cor.cc.sn[non.na,non.na], col = cols, zlim = c(-1,1),main = 'CC vs SN')
}

#MOV vs all
quartz(width = 10, height = 10)
par(mfcol=c(length(sessions),3), mar = c(4,4,2,4))
for(i in 1:length(sessions)){
  non.na <- which(rowSums(!is.na(cor.mats[[i]]$cor.cc.cc)) >0)
  image.plot(cor.mats[[i]]$cor.mov.mov[non.na,non.na], col = cols, zlim = c(-1,1), main = 'MOV vs MOV')
  image.plot(cor.mats[[i]]$cor.cc.mov[non.na,non.na], col = cols, zlim = c(-1,1), main = 'CC vs MOV')
  image.plot(cor.mats[[i]]$cor.mov.sn[non.na,non.na], col = cols, zlim = c(-1,1),main = 'MOV vs SN')
}

#SN vs all
quartz(width = 10, height = 10)
par(mfcol=c(length(sessions),3), mar = c(4,4,2,4))
for(i in 1:length(sessions)){
  non.na <- which(rowSums(!is.na(cor.mats[[i]]$cor.cc.cc)) >0)
  image.plot(cor.mats[[i]]$cor.sn.sn[non.na,non.na], col = cols, zlim = c(-1,1), main = 'SN vs SN')
  image.plot(cor.mats[[i]]$cor.cc.sn[non.na,non.na], col = cols, zlim = c(-1,1), main = 'CC vs SN')
  image.plot(cor.mats[[i]]$cor.mov.sn[non.na,non.na], col = cols, zlim = c(-1,1),main = 'MOV vs SN')
}

#  How does call rate vary with nearest neighbor distance and nearest neighbor call rate?

#get labels
cc.rate.labs <- paste0(cc.rate.bins[1:length(cc.rate.bins)-1])
cc.rate.labs[length(cc.rate.labs)] <- paste0(cc.rate.labs[length(cc.rate.labs)],'+')
move.rate.labs <- paste0(move.rate.bins[1:length(move.rate.bins)-1])
move.rate.labs[length(move.rate.labs)] <- paste0(move.rate.labs[length(move.rate.labs)],'+')
sn.rate.labs <- paste0(sn.rate.bins[1:length(sn.rate.bins)-1])
sn.rate.labs[length(sn.rate.labs)] <- paste0(sn.rate.labs[length(sn.rate.labs)],'+')
nn.dist.labs <- paste0(nn.dist.bins[1:(length(nn.dist.bins)-1)])
nn.dist.labs[length(nn.dist.labs)] <- paste0(nn.dist.labs[length(nn.dist.labs)],'+')

quartz(width = 16, height = 12)
par(mfcol=c(length(sessions), 3), cex.lab =1.5, cex.axis = 1.5, cex.main=2, mar = c(5,5,5,5))
for(sess in 1:length(sessions)){
  
  #cc 
  image.plot(call.rate.vs.nn.mats[[sess]]$cc.rate.vs.nn.mat, col = viridis(256), xaxt = 'n', yaxt = 'n', xlab = 'NN call rate', ylab = 'NN dist', main = paste(sessions[sess],'CC',sep=' - '))
  axis(side = 1, at = seq(0,1,length.out=length(cc.rate.labs)),labels = cc.rate.labs)
  axis(side = 2, at = seq(0,1,length.out=length(nn.dist.labs)),labels = nn.dist.labs)
  
  #move
  image.plot(call.rate.vs.nn.mats[[sess]]$move.rate.vs.nn.mat, col = viridis(256), xaxt = 'n', yaxt = 'n', xlab = 'NN call rate', ylab = 'NN dist', main = paste(sessions[sess],'MOV',sep=' - '))
  axis(side = 1, at = seq(0,1,length.out=length(move.rate.labs)),labels = move.rate.labs)
  axis(side = 2, at = seq(0,1,length.out=length(nn.dist.labs)),labels = nn.dist.labs)
  
  #sn
  image.plot(call.rate.vs.nn.mats[[sess]]$sn.rate.vs.nn.mat, col = viridis(256), xaxt = 'n', yaxt = 'n', xlab = 'NN call rate', ylab = 'NN dist', main = paste(sessions[sess],'SN',sep=' - '))
  axis(side = 1, at = seq(0,1,length.out=length(sn.rate.labs)),labels = sn.rate.labs)
  axis(side = 2, at = seq(0,1,length.out=length(nn.dist.labs)),labels = nn.dist.labs)
  
}
