#This script uses cross-correlation functions to explore call/response dynamics as a function time and space

#It works as follows:
#Make a table of all calls within a given call category (set of calls), identify the one who gave the call as the 'caller'
#For each call, identify all other individuals ('responders') who were present (calls labeled) on that day, note their distance from the caller
#Also identify the distance between the caller and responder at the time when the call was given by the caller
#Save this table as callresp

#Then for each of row in this table, look at the call sequence of the 'responder' over time, lined up so that t = 0 = time of the 'caller' call
#Run a kernel over this time series to smooth it with a bandwidth bw = .1 sec
#Save the output to a matrix as callresp.seqs, such that the indices in the table match the indices in the matrix

#Finally, get the mean across all (or a subset) of rows of the matrix for different conditions:
# Different distance bins (aggregate only data from caller-responder pairs within a range of distances from one another)
# Different individuals + distance bins (aggregate only data from caller-responder pairs where the caller was a certain individual, and distance < 3 m since this was discovered to be the relevant range from the distance analysis)

#----------------PARAMETERS--------------------

#whether to save output
save.output <- T

#directory where data is stored
datadir <- '~/Dropbox/meerkats/meerkats_shared/data' 

#directory where code is stored for this project
codedir <- '~/Dropbox/code_ari/meerkat_vocal_interactions'

#directory of where to save results
savedir <- '~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/data/call_response'

#filename of general meerkat functions
general.funcs.filename <- '~/Dropbox/code_ari/move_comm_analysis/general/meerkat_functions.R'

#bandwidth of smoothing kernel (default 0.1)
bw <- .05

#maximum time lag to consider (time since another individual called)
max.lag <- 30 

#time step to use for the sequence of times
step <- .02

#list of call types to include in the set of calls by the initial caller (which determines the 0 point of the correlogram)
caller.calltypes <- c('cc') 

#list of call types of include in the set of calls by the responder (determines the curve in the correlogram)
responder.calltypes <- c('cc')

#what to trigger on (call begin or call end)
trigger.on <- 'begin'

#store parameters in a named list
params <- list()
params$datadir <- datadir
params$codedir <- codedir
params$bw <- bw
params$max.lag <- max.lag
params$step <- step
params$caller.calltypes <- caller.calltypes
params$responder.calltypes <- responder.calltypes
params$trigger.on <- trigger.on

#whether to instead perform an analysis of whether an individual repeats its call in a certain time window after its initial call
#as a function of whether there have been any calls from nearby individuals (within repeat.dist.thresh meters) in the intervening time
#all sequences after a call are removed if there has been a call by anyone in the past repeat.time.thresh seconds
repeat.self.analysis <- FALSE
repeat.time.thresh <- 10
repeat.dist.thresh <- 10

#---------------------SETUP-----------------------

#set time zone to UTC
Sys.setenv(TZ='UTC')

#list of group years
groupyears <- c('HM2017', 'HM2019', 'L2019')

#file names
audio.file <- 'full_labelfile_conflicts_resolved.csv'
gps.files <- paste(groupyears, 'COORDINATES_all_sessions.RData', sep = '_')

if(save.output){
  savename <- paste0('callresp_', caller.calltypes, '_', responder.calltypes, '_bw', bw, '.RData')
}

#------------------LIBRARIES----------------------

#libraries
library(jcolors)
library(fitdistrplus)
library(viridis)
library(fields)

#------------------LOAD DATA-----------------------

timestamp()
print('loading data')

#set working directory
setwd(datadir)

#load audio data
calls.all <- read.csv(audio.file, header=T, sep='\t', stringsAsFactors=F)

#load gps data
for(i in 1:length(gps.files)){
  load(gps.files[i])
}

#------------------FUNCTIONS------------------------
#source functions
setwd(codedir)

source(general.funcs.filename)

#--------------------PREPROCESS-------------------------------

timestamp()
print('preprocessing')

#check parameters
if(!(trigger.on %in% c('begin','end'))){
  stop('trigger.on parameter should be either begin or end')
}

#convert times to POSIXlt
calls.all$t0 <- as.POSIXlt(calls.all$t0GPS_UTC, tz = 'UTC')
calls.all$tf <- as.POSIXlt(calls.all$tEndGPS_UTC, tz = 'UTC')

#remove double-marked calls (same start and end time and same individual calling)
dups <- which(duplicated(cbind(calls.all$ind,as.numeric(calls.all$t0),as.numeric(calls.all$tf))))
if(length(dups)>0){
  calls.all <- calls.all[-dups,]
}

#get unique dates
dates <- unique(calls.all$date)

#indicate which calls are caller call type list, and which are in repsonder call type list
calls.all$isCallerCallType <- F
calls.all$isResponderCallType <- F
calls.all$isCallerCallType[which(calls.all$callType %in% caller.calltypes)] <- T
calls.all$isResponderCallType[which(calls.all$callType %in% responder.calltypes)] <- T

#---------------------------MAIN-------------------------------
#time sequence bins
tseq <- seq(-max.lag, max.lag, step)

#get together a table of all calls we are using in the analysis t0 and tf times from all individuals, as well as date and the start and end time of the labeled sequence
calls.use <- data.frame()

for(d in 1:length(dates)){
  date <- dates[d]
  
  #get calls for that date
  calls.date <- calls.all[which(calls.all$date == date),]
  
  #get start and stop times for times when all are labeled - need to fix this a bit because sometimes there are multiiple start markers (e.g. VHMF010 on 20190712)
  starts <- calls.date[which(calls.date$callType %in% c('start','START')),c('ind','t0')]
  ends <- calls.date[which(calls.date$callType %in% c('stop','end','STOP','END')),c('ind','t0')]
  
  #for now just take the minimum 'start' for each individual as its start marker and the max end marker as its end time
  starts <- aggregate(starts$t0,by = list(starts$ind), min)
  ends <- aggregate(ends$t0,by = list(ends$ind), max)
  
  #find the latest start and earliest end times to be the group start and end time
  start.all <- max(starts$x)
  end.all <- min(ends$x)
  
  #make a simple table that just contains the cc's from each individual at that date and their t0 and tf times (as numeric)
  calls.use.date <- calls.date[which((calls.date$isCallerCallType | calls.date$isResponderCallType) & calls.date$pred_focalType == 'F'),]
  
  #append to the larger table
  if(nrow(calls.use.date)>0){
    calls.use.date <- calls.use.date[,c('date','ind','callType','isCallerCallType','isResponderCallType','t0GPS_UTC','tEndGPS_UTC','t0','tf')]
    calls.use.date$t0 <- as.numeric(calls.use.date$t0)
    calls.use.date$tf <- as.numeric(calls.use.date$tf)
    calls.use.date$start <- as.numeric(start.all)
    calls.use.date$end <- as.numeric(end.all)
    
    calls.use <- rbind(calls.use, calls.use.date)
  }
}

colnames(calls.use) <- c('date','caller','callType','isCallerCallType','isResponderCallType','t0GPS','tfGPS','t0','tf','start','end')

#create a list of which individuals were labeled on which date
inds.present <- list()
for(i in 1:length(dates)){
  inds.present[[i]] <- unique(calls.all$ind[which(calls.all$date == dates[i])])
}

#create an expanded data frame that contains rows for each other 'responder' individual (not the caller)
timestamp()
print('generating call-response table')
callresp <- data.frame()
for(i in 1:nrow(calls.use)){
  
  if(i %% 10000 == 0){
    print(paste(i,'/',nrow(calls.use)))
  }
  
  row <- calls.use[i,]
  
  if(row$isCallerCallType == T){
    
    date.idx <- which(dates == calls.use$date[i])
    inds <- inds.present[[date.idx]]
    
    new.rows <- row[rep(1,length(inds)),]
    new.rows$responder <- inds
    
    callresp <- rbind(callresp, new.rows)
  }
  
}

#add a column with the correspoding groupyear to the callresp table
groupyear.dates <- list()
for(i in 1:length(groupyears)){
  gy <- groupyears[i]
  gy_tl <- eval(as.name(paste(gy,'timeLine', sep = '_')))
  gy_dates <- unique(date(gy_tl))
  gy_date_char <- gsub('-','', as.character(gy_dates))
  groupyear.dates[[i]] <- gy_date_char
}

#add a column to the callresp table with the groupyear (to know which GPS data to use)
callresp$groupyear <- NA
for(i in 1:length(groupyears)){
  idxs <- which(callresp$date %in% groupyear.dates[[i]])
  callresp$groupyear[idxs] <- groupyears[i]
}

#get distance between caller and responder at time of the call (t0)
#loop over group years
timestamp()
print('computing distances between individuals across all recording periods')
#TODO: First get time indexes, only do once to speed up code
callresp$distance <- NA
for(g in 1:length(groupyears)){
  
  print(groupyears[g])
  
  #get indices associated with that group year, and associated timeline
  idxs <- which(callresp$groupyear == groupyears[g])
  timeLine <- eval(as.name(paste(groupyears[g],'timeLine', sep = '_')))
  timeline.numeric <- as.numeric(as.POSIXct(timeLine, tz = 'UTC'))
  
  allX <- eval(as.name(paste(groupyears[g], 'allX', sep = '_')))
  allY <- eval(as.name(paste(groupyears[g], 'allY', sep = '_')))
  indInfo <- eval(as.name(paste(groupyears[g],'indInfo', sep = '_')))
  
  #loop over all indices associated with that group year and get response sequences
  for(i in idxs){
    caller <- callresp$caller[i]
    responder <- callresp$responder[i]
    
    caller.idx <- match(caller,indInfo$code)
    responder.idx <- match(responder,indInfo$code)
    
    t0 <- round(callresp$t0[i]) ##
    t0.idx <- which(timeline.numeric==t0)
    
    if(length(t0.idx)!=0){
      callresp$distance[i] <- sqrt((allX[caller.idx,t0.idx] - allX[responder.idx,t0.idx])^2 + (allY[caller.idx,t0.idx] - allY[responder.idx,t0.idx])^2)
    }
    
  }
  
}

#create a matrix to hold call/response sequences
callresp.seqs <- matrix(data=NA,nrow=nrow(callresp),ncol = length(tseq))

#loop through and get individual call patterns starting at tf
timestamp()
print('gathering vocal sequences across all recording periods')
for(i in 1:nrow(callresp.seqs)){
  
  if(trigger.on == 'begin'){
    zero.time <- callresp$t0[i]
  } else{
    zero.time <- callresp$tf[i]
  }
  responder <- callresp$responder[i]
  caller <- callresp$caller[i]
  
  #if too close to end, don't use
  if(((zero.time - max.lag) > callresp$start[i]) & ((zero.time + max.lag) < callresp$end[i])){
    
    #get call times of the responder (only of the correct call type, and removing the current call)
    current.call.idx <- which(calls.use$caller == caller & calls.use$t0 == callresp$t0[i] & calls.use$tf == callresp$tf[i])
    if(length(current.call.idx)==0){
      stop('could not find current call in the table')
    } 
    if(length(current.call.idx) > 1){
      warning('found two calls with same label, start, and stop time')
    }
    calls.use.removeocurrent <- calls.use[-current.call.idx,]
    call.times.foc <- calls.use.removeocurrent$t0[which(calls.use.removeocurrent$isResponderCallType & calls.use.removeocurrent$caller == responder)]
    
    dt.foc <- call.times.foc - zero.time
    dt.foc <- dt.foc[which(abs(dt.foc) < max.lag)] #get rid of way early or way late calls to reduce compute time
    call.seq <- rep(0, length(tseq))
    
    if(length(dt.foc)>0){
      for(j in 1:length(dt.foc)){
        
        #for self-response version, don't include the current call
        if((self.responses | repeat.self.analysis) & dt.foc[j] == 0){
          call.seq <- call.seq
        } else{
          call.seq <- call.seq + dnorm(tseq, mean = dt.foc[j], sd = bw)
        }
        
        
      }
    }
    
    callresp.seqs[i,] <- call.seq
    
  }
}

timestamp()

save(list = c('callresp','callresp.seqs'))

print('plotting results')

#if performing an analysis of whether you repeat a call depending on whether others have replied, need to compute some extra things
#TODO: CHECK THIS AFTER UPDATES 2021-05-07
if(repeat.self.analysis == T){
  
  #----SETUP---
  #first the time of the previous call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
  callresp.filt$prevcall.t0 <- NA
  callresp.filt$prevcall.caller <- NA
  callresp.filt$prevcall.type <- NA
  callresp.filt$n.inds.in.range <- NA
  
  #time of the next call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
  callresp.filt$nextcall.t0 <- NA
  callresp.filt$nextcall.caller <- NA
  callresp.filt$nextcall.type <- NA
  
  for(i in 1:nrow(callresp.filt)){
    
    #get time of the focal call
    t0 <- callresp.filt$t0[i] 
    
    #time index in timeLine
    t0.round <- round(t0)
    t.idx <- which(timeline.numeric == t0.round)
    
    #get individual index
    ind.idx <- which(indInfo$code == callresp.filt$caller[i])
    
    #distances to all other individuals present
    dists <- sqrt((allX[,t.idx] - allX[ind.idx,t.idx])^2 + (allY[,t.idx] - allY[ind.idx,t.idx])^2)
    
    #distance to self is NA
    dists[ind.idx] <- NA
    
    #get the individuals that are in range (wihtin repeat.dist.thresh)
    inds.in.range.idxs <- which(dists < repeat.dist.thresh)
    
    #store number of individuals in range
    callresp.filt$n.inds.in.range[i] <- length(inds.in.range.idxs)
    
    #if there are some indivudals in range, find the most recent call from any of them of any type, and next call of a type in responder.calltypes  
    if(length(inds.in.range.idxs)>0){
      
      #get individuals in range at the time of the focal call
      inds.in.range <- indInfo$code[inds.in.range.idxs]
      
      #get all calls from the individuals in range
      calls.inds.in.range <- calls.all[which(calls.all$ind %in% inds.in.range & calls.all$isCall==1 & calls.all$pred_focalType=='F'),]
      
      #get most recent previous call from the individuals in range (of any type)
      calls.before <- calls.inds.in.range$t0[which(calls.inds.in.range$t0 < t0)]
      if(length(calls.before) > 0){
        prev.call.time.inds.in.range <- max(calls.before)
        callresp$prevcall.t0[i] <- as.numeric(prev.call.time.inds.in.range, tz = 'UTC')
        prev.call.inds.in.range.idx <- which(calls.inds.in.range$t0 == prev.call.time.inds.in.range)[1]
        callresp$prevcall.caller[i] <- calls.inds.in.range$ind[prev.call.inds.in.range.idx]
        callresp$prevcall.type[i] <- calls.inds.in.range$callSimple[prev.call.inds.in.range.idx]
      }
      
      #get next call from individuals in range (of types in responder.calltypes)
      correct.calls.after <- calls.inds.in.range$t0[which(calls.inds.in.range$t0 > t0 & calls.inds.in.range$callType %in% responder.calltypes)]
      if(length(correct.calls.after)>0){
        correct.calls.after.numeric <- as.numeric(correct.calls.after, tz = 'UTC')
        next.call.time.inds.in.range <- min(correct.calls.after.numeric)
        callresp.filt$nextcall.t0[i] <- as.numeric(next.call.time.inds.in.range)
        next.call.inds.in.range.idx <- which(calls.inds.in.range$t0 == next.call.time.inds.in.range)[1]
        callresp.filt$nextcall.caller[i] <- calls.inds.in.range$ind[next.call.inds.in.range.idx]
        callresp.filt$nextcall.type[i] <- calls.inds.in.range$callSimple[next.call.inds.in.range.idx]
      }
    }
  }
  
  callresp.filt$dt.prev <- callresp.filt$t0 - callresp.filt$prevcall.t0
  callresp.filt$dt.next <- callresp.filt$nextcall.t0 - callresp.filt$t0
}

if(save.output){
  save(file = paste0(savedir, '/', savename), list = c('callresp','callresp.seqs','tseq', 'params'))
}

#---------------------PLOTTING-------------------

#ALL INDIVIDUALS - DIFFERENT DISTANCE RANGES

#Make a plot of cross-correlogram across all individuals

#collect the data
mean.call.rates <- matrix(NA,nrow=length(dist.bins)-1, ncol = length(tseq))
for(i in 2:length(dist.bins)){
  idxs <- which(callresp.filt$distance >= dist.bins[i-1] & callresp.filt$distance < dist.bins[i])
  mean.call.rates[i-1,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
}

#make the plot
quartz(height = 8, width = 12)
par(mfrow=c(1,length(dist.bins)-1))
par(mar=c(8,6,3,1))
ymax <- max(mean.call.rates)+.01
cols <- viridis(nrow(mean.call.rates))
for(i in 1:nrow(mean.call.rates)){
  plot(NULL, xlim=c(-2,2),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5, main = paste(dist.bins[i],'-', dist.bins[i+1],'m',sep=' '))
  abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
  abline(v=0, lty=1, col = 'black')
  lines(tseq, mean.call.rates[i,], col = cols[i], lwd = 2, type = 'l', )
}

#make a "spiking neurons" plot
i <- 2
n.neurons <- 200
xmax <- 10
isna <- is.na(rowSums(callresp.seqs))
idxs <- which((callresp.filt$distance >= dist.bins[i-1]) & (callresp.filt$distance < dist.bins[i]) & !isna)
idxs.sub <- sample(idxs, n.neurons)

quartz(height = 8, width = 4)
plot(NULL, xlim = c(-xmax, xmax), ylim = c(0, n.neurons), xlab = "Time (sec)", ylab = 'Calls')
for(j in 1:n.neurons){
  lines(tseq, callresp.seqs[idxs.sub[j],] + j )
}
abline(v = 0, col = 'blue')



# #SUBSET BY CALLER
# 
# #subset by caller (here only use 0-3 meter length scale as this is what was discovered before)
# mean.call.rates <- matrix(NA,nrow=length(inds), ncol = length(tseq))
# n.samps <- rep(NA, length(inds))
# for(i in 1:length(inds)){
#   idxs <- which(callresp.filt$caller==inds[i] & callresp.filt$distance < 3)
#   mean.call.rates[i,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
#   n.samps[i] <- length(idxs)
# }
# 
# quartz(height=12,width = 12)
# par(mfrow=c(ceiling(length(inds)/3),3))
# par(mar=c(6,5,1,1))
# 
# #set up color palette
# cols <- jcolors('pal8')
# ymax <- max(mean.call.rates,na.rm=T)+.01
# for(i in 1:length(inds)){
#   plot(tseq, mean.call.rates[i,], col = cols[i], lwd = 2,type='l',xlim=c(-5,5),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5)
#   abline(v=0, lty=2)
#   legend(x = 'topleft',legend=inds[i],col=cols[i],lwd=1.5, cex=1.5)
# }
# 
# #if performing analysis of whether calls are repeated when a reply is vs isn't given, make this plot
# #TODO: Generalize to include the new data
# if(repeat.self.analysis){
#   ##bins to measure individual call rate over (after a focal call), get response sequences only after t0 = 0
#   bins <- tseq[which(tseq >= 0)]
#   
#   #get only positive time lags
#   callresp.seqs.repeat <- callresp.seqs[,which(tseq >= 0)]
#   
#   #loop over bins to compute call rates when there has or hasn't been a response
#   callrate.resp <- callrate.nonresp <- callrate.nonbrs <- n.resp <- n.nonresp <- n.nonbr <- rep(NA, length(bins)-1)
#   for(i in 2:length(bins)){
#     response.idxs <- which(callresp$dt.next < bins[i] & callresp$dt.prev > repeat.time.thresh)
#     nonresponse.idxs <- which(callresp$dt.next > bins[i] & callresp$dt.prev > repeat.time.thresh)
#     noneighbor.idxs <- which(callresp$n.inds.in.range == 0)
#     
#     callrate.resp[i] <- mean(callresp.seqs.repeat[response.idxs,i],na.rm=T)
#     callrate.nonresp[i] <- mean(callresp.seqs.repeat[nonresponse.idxs,i],na.rm=T)
#     callrate.nonbrs[i] <- mean(callresp.seqs.repeat[noneighbor.idxs,i],na.rm=T)
#     n.resp[i] <- length(response.idxs)
#     n.nonresp[i] <- length(nonresponse.idxs)
#     n.nonbr[i] <- length(noneighbor.idxs)
#   }
#   
#   quartz(height = 12, width = 8)
#   par(mfrow=c(2,1))
#   plot(bins, callrate.resp, type='l', col = 'red', lwd = 2, xlab = 'Time after call (sec)', ylab = c('Call rate'), cex.lab = 1.5, main = paste('Min silence = ', repeat.time.thresh, ' | Max dist = ', repeat.dist.thresh), ylim = c(0,0.2))
#   lines(bins, callrate.nonresp, col = 'blue', lwd = 2, lty = 1)
#   lines(bins, callrate.nonbrs, col = 'black', lwd = 2, lty = 1)
#   legend('topright', legend = c('reply', 'no reply','no neighbor'), col = c('red','blue','black'), lty = c(1,1), cex = 1.5)
#   plot(bins, n.resp, type = 'l', lwd = 2, xlab = 'Time (sec)', ylab = 'Number of events', col = 'red', cex = 2, cex.lab = 1.5, ylim = c(0, max(n.resp + n.nonresp + n.nonbr,na.rm=T)))
#   lines(bins,n.nonresp, lwd = 2, col = 'blue')
#   lines(bins,n.nonbr, lwd = 2, col = 'black')
# }