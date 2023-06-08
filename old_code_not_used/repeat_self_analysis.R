
#------PARAMETERS------
repeat.dist.thresh <- 5
n.boots <- 100

#-----FILENAME------
filename <- '~/Dropbox/meerkats/results/call_interactions/callresp_cc_cc_bw0.05_30sec.RData'

#-----LOAD------
load(filename)

#----SETUP---

#remove callresp.seqs first half (the sequence from before a given call) 
#also reduce the temporal resolution to save compute time
nT <- length(tseq)
res_idxs <- seq(from=((nT-1)/2+1),to=nT, by = 5)
callresp.seqs <- callresp.seqs[,res_idxs]
tseq <- tseq[res_idxs]

#break into self-responses and other responses
self_idxs <- which(callresp$caller==callresp$responder)
other_idxs <- which(callresp$caller!=callresp$responder)
callresp_self <- callresp[self_idxs,]
callresp_other <- callresp[other_idxs,]
callresp.seqs_self <- callresp.seqs[self_idxs,]
callresp.seqs_other <- callresp.seqs[other_idxs,]
rm(list=c('callresp','callresp.seqs'))

#retain only individuals within 5 m
callresp_other_close <- callresp_other[which(callresp_other$distance <= repeat.dist.thresh),]

#find how many individuals were within the distance threshold
#also find the time of the first response by an individual within the distance threshold (if any)
callresp_self$t_first_response <- NA
callresp_self$inds_in_range <- NA
rows <- seq(1, nrow(callresp_self))

timestamp()
for(i in rows){
  
  if(i %% 1000 == 1){
    print(paste0(i,'/',length(rows)))
  }
  
  #time of the initial call
  t0 <- callresp_self$t0[i]
  
  #caller of the initial call
  caller <- callresp_self$caller[i]
  
  #how many individuals were located within 5 m of the caller at the time of the original call?
  responder_idxs <- which(callresp_other_close$caller == caller & callresp_other_close$t0 == t0)
  n_in_range <- length(responder_idxs)
  callresp_self$inds_in_range[i] <- n_in_range
  
  #get all calls associated with other callers which are after t0 and within 5 m of the original caller
  later_call_idxs <- which(callresp_other_close$t0 >= t0 & callresp_other_close$caller != caller)
  
  #find the first call from amongst those within 5 m of the original caller (if no calls, return NA)
  if(length(later_call_idxs)==0){
    next
  }
  first_call_time <- min(callresp_other_close$t0[later_call_idxs])
  dt <- first_call_time-t0
  callresp_self$t_first_response[i] <- dt
  
}
timestamp()

#get index (in tseq vector) after which the individual has been responded to
#this gives warnings - it's OK
callresp_self$tidx_first_response <- sapply(callresp_self$t_first_response, FUN = function(x){return(min(which(tseq>=x)))})
callresp_self$tidx_first_response[which(is.infinite(callresp_self$tidx_first_response))] <- NA

#create callresp.seqs_self for when an individual does and does not hear a response within 5 m
#first just copy over the callresp.seqs_self matrix to both matrices
callseq_beforeresp <- callseq_afterresp <- callresp.seqs_self

#then fill with NAs AFTER the first response for callseq_beforeresp and BEFORE the first response for callseq_afterresp
maxT <- max(tseq)
max_tidx <- length(tseq)
for(i in rows){
  
  #if there was never a response, replace all values in that row for both matrices with NAs
  if(is.na(callresp_self$t_first_response[i])){
    callseq_afterresp[i,] <- NA
    callseq_beforeresp[i,] <- NA
    next
  }
  
  #if the first response came after maxT, consider everything to be pre-response
  if(callresp_self$t_first_response[i] >= maxT){
    callseq_afterresp[i,] <- NA
    next
  }
  
  #get the index of the tiem of first response
  tidx <- callresp_self$tidx_first_response[i]
  
  #otherwise, fill in all values before tidx with NAs for callseq_afterresp
  #and fill in all values after tidx with NAs for callseq_beforeresp
  callseq_beforeresp[i,tidx:max_tidx] <- NA
  callseq_afterresp[i,1:(tidx-1)] <- NA
  
}

#either use all data or make some subset
idxs_to_remove <- which(callresp_self$inds_in_range==0) #only when others are around
idxs_to_remove <- NULL #don't remove anything

#put NAs in the idxs we decided to remove
if(length(idxs_to_remove)>0){
  callseq_beforeresp[idxs_to_remove,] <- NA
  callseq_afterresp[idxs_to_remove,] <- NA
}

#get means

means_before <- colMeans(callseq_beforeresp,na.rm=T)
means_after <- colMeans(callseq_afterresp,na.rm=T)

#bootstrap by groupyear and date
callresp_self$groupyeardate <- paste0(callresp_self$groupyear,'_',callresp_self$date)
groupyeardates <- unique(callresp_self$groupyeardate)
means_before_boot <- means_after_boot <- matrix(nrow=length(means_before),ncol=n.boots)
for(i in 1:n.boots){
  print(i)
  groupyeardates_boot <- sample(groupyeardates,replace=T)
  idxs <- c()
  for(j in 1:length(groupyeardates_boot)){
    idxs <- c(idxs, which(callresp_self$groupyeardate == groupyeardates_boot[j]))
  }
  means_before_boot[,i] <- colMeans(callseq_beforeresp[idxs,],na.rm=T)
  means_after_boot[,i] <- colMeans(callseq_afterresp[idxs,],na.rm=T)
  
}

#plot
maxy <- max(max(means_before_boot,na.rm=T),max(means_after_boot,na.rm=T))
quartz()
plot(NULL, xlim=c(0,20),ylim=c(0,maxy),xlab = 'Time after initial call (sec)',ylab = 'Self-reply rate')

#for(i in 1:n.boots){
#  lines(tseq,means_before_boot[,i],col='red',lwd=.2)
#  lines(tseq,means_after_boot[,i],col='blue',lwd=.2)
#}

uppers_before <- apply(means_before_boot, 1, function(x){return(quantile(x,0.975,na.rm=T))})
uppers_after <- apply(means_after_boot, 1, function(x){return(quantile(x,0.975,na.rm=T))})
lowers_before <- apply(means_before_boot, 1, function(x){return(quantile(x,0.025,na.rm=T))})
lowers_after <- apply(means_after_boot, 1, function(x){return(quantile(x,0.025,na.rm=T))})

polygon(c(tseq,rev(tseq)), c(uppers_before, rev(lowers_before)), border=NA, col = alpha('red',0.2))
polygon(c(tseq,rev(tseq)), c(uppers_after, rev(lowers_after)), border=NA, col = alpha('blue',0.2))

lines(tseq,means_before,col='red',lwd=3)
lines(tseq,means_after,col='blue',lwd=3)
legend('bottomright', legend = c('Response within 5 m', 'No response within 5 m'), col = c('blue','red'),lwd=c(3,3))

