#This code does not currently work!
#It was taken from the original vocal_interactions_across_space_and_time script and needs to be integrated into a subsequent script to work (after computing call response sequences)

#------PARAMETERS------
repeat.dist.thresh <- 10

#-----FILENAME------
filename <- '~/Dropbox/meerkats/meerkats_shared/ari/vocal_interactions/data/call_response/callresp_cc_cc_bw0.05.RData'

#-----LOAD------
load(filename)

#----SETUP---
#first the time of the previous call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
callresp$prevcall.t0 <- NA
callresp$prevcall.caller <- NA
callresp$prevcall.type <- NA
callresp$n.inds.in.range <- NA

#time of the next call (by anyone within repeat.dist.thresh of the caller), the caller id, and the call type
callresp$nextcall.t0 <- NA
callresp$nextcall.caller <- NA
callresp$nextcall.type <- NA

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

