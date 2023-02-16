library(viridis)
library(scales)

#-------------------FILENAME--------------------
filename <- '~/Dropbox/meerkats/results/call_interactions/callresp_cc_cc_bw0.1.RData'

#data directory where ind info is stored
#ind_info_dir <- '/Volumes/EAS_shared/meerkat/working/METADATA/'
ind_info_dir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/METADATA/'

#------------------PARAMETERS------------------
dist.bins <- c(0,2,5,10,50)
#dist.bins <- c(0,.1,.5,1,2,3,5,10,50)
n.boots <- 100 #number of bootstraps to do for error bars
gy <- NULL #which group year to use (if NULL, use all groups together)

#plots to do
plot.spiking.neurons <- F
plot.call.resp.all <- T
plot.self.resp.all <- F
plot.by.caller <- F

#caller age classes
adult_classes <- c('DominantF','DominantM','Yearling','Sub-Adult','Adult')
juv_classes <- c('Juvenile')

#--------------------HELPER FUNCS---------------
#Calculate mean call rates from data frame of call response data and time sequences
#INPUTS:
# callresp: data frame containing info on call responses (from get_call_response_sequences script)
# callresp.seqs: matrix containing the actual time series data (from get_call_response_sequences script)
# gy: groupyear, currently 'HM2017','HM2019', and 'L2019' are supported. If NULL, use all groups together
# tseq: vector sequence of times (from get_call_response_sequences script)
# dist.bins: vector of distance bins (default to 0,2,5,10,50)
# adult_classes: age classes to use (default to c('DominantF','DominantM','Yearling','Sub-Adult','Adult'))
get_mean_call_rates <- function(callresp, callresp.seqs, gy, tseq, 
                                dist.bins = c(0,2,5,10,50), 
                                adult_classes = c('DominantF','DominantM','Yearling','Sub-Adult','Adult')){
  mean.call.rates <- matrix(NA,nrow=length(dist.bins)-1, ncol = length(tseq))
  for(i in 2:length(dist.bins)){
    if(!is.null(gy)){
      idxs <- which((callresp$distance >= dist.bins[i-1]) &
                    (callresp$distance < dist.bins[i]) & 
                    (callresp$caller != callresp$responder) & 
                    (callresp$age_caller %in% adult_classes) & 
                    (callresp$age_responder %in% adult_classes)&
                    callresp$groupyear == gy)
    } else{
      idxs <- which((callresp$distance >= dist.bins[i-1]) &
                      (callresp$distance < dist.bins[i]) & 
                      (callresp$caller != callresp$responder) & 
                      (callresp$age_caller %in% adult_classes) & 
                      (callresp$age_responder %in% adult_classes))
      }
    mean.call.rates[i-1,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
    }
  return(mean.call.rates)
}

#--------------------LOAD DATA------------------
load(filename)

setwd(ind_info_dir)
#individual info
groupyears <- unique(callresp$groupyear)
ind_info <- data.frame()
for(g in 1:length(groupyears)){
  groupyear <- groupyears[g]
  tmp <- read.csv(paste(groupyear,'_INDIVIDUAL_INFO.txt',sep=''),sep='\t')
  tmp$groupyear <- groupyear
  tmp$groupyear_ind <- paste(tmp$groupyear, tmp$code, sep='_')
  ind_info <- rbind(ind_info, tmp)
}

#add columns to specify age class
callresp$groupyear_caller <- paste(callresp$groupyear,callresp$caller, sep='_')
callresp$groupyear_responder <- paste(callresp$groupyear, callresp$responder, sep = '_')
callresp$age_caller <- ind_info$status[match(callresp$groupyear_caller, ind_info$groupyear_ind)]
callresp$age_responder <- ind_info$status[match(callresp$groupyear_responder, ind_info$groupyear_ind)]

#---------------------PLOTTING-------------------

#---PLOT 1: Spiking neurons example
if(plot.spiking.neurons){
  i <- 2
  n.neurons <- 10
  xmax <- 10
  isna <- is.na(rowSums(callresp.seqs))
  idxs <- which((callresp$distance >= dist.bins[i-1]) & (callresp$distance < dist.bins[i]) & !isna)
  idxs.sub <- sample(idxs, n.neurons)
  
  quartz(height = 8, width = 4)
  plot(NULL, xlim = c(-xmax, xmax), ylim = c(0, n.neurons), xlab = "Time (sec)", ylab = 'Calls')
  for(meerj in 1:n.neurons){
    lines(tseq, callresp.seqs[idxs.sub[j],] + j )
  }
  abline(v = 0, col = 'blue')
}

#CALL RESPONSE PLOTS
if(plot.call.resp.all){
  
  #---PLOT 2: Call response dynamics on average, between two adult individuals
  #At different spatial scales
  #collect the data

  #get mean call rates for each distance bin and time bin
  mean.call.rates <- get_mean_call_rates(callresp, callresp.seqs, gy, tseq, dist.bins, adult_classes)
  
  #get mean call rates in bootstrapped data
  #matrix to hold output
  mean.call.rates.boot <- array(NA,dim = c(length(dist.bins)-1,length(tseq), n.boots))
  
  #unique identifier for selecting relevant sequences in bootstrapping
  callresp$groupyeardate <- paste0(callresp$groupyear,'_',callresp$date)
  groupyeardates <- unique(callresp$groupyeardate)
  
  #bootstrap by selecting same number of days as we have in real dataset with replacement
  for(b in 1:n.boots){
    print(paste0(b,'/',n.boots))
    timestamp()
    
    #bootstrap for error bars - either by call or by day/group
    groupyeardates.boot <- sample(groupyeardates, replace=T)
    groupyeardates.boot.idxs <- c()
    for(bb in 1:length(groupyeardates.boot)){
      groupyeardates.boot.idxs <- c(groupyeardates.boot.idxs, which(callresp$groupyeardate == groupyeardates.boot[bb]))
    }
    callresptmp <- callresp[groupyeardates.boot.idxs,]
    callrespseqstmp <- callresp.seqs[groupyeardates.boot.idxs,]
    mean.call.rates.boot[,,b] <- get_mean_call_rates(callresptmp, callrespseqstmp, gy, tseq, dist.bins, adult_classes)
  }
  
  #get upper and lowers (95% interval)
  uppers.boot <- apply(mean.call.rates.boot, c(1,2), FUN = function(x){return(quantile(x,0.975,na.rm=T))})
  lowers.boot <- apply(mean.call.rates.boot, c(1,2), FUN = function(x){return(quantile(x,0.025,na.rm=T))})
  
  #make the plot
  quartz(height = 6, width = 12)
  par(mfrow=c(1,length(dist.bins)-1))
  par(mar=c(6,6,3,1))
  ymax <- max(mean.call.rates)*1.3
  cols <- viridis(nrow(mean.call.rates))
  for(i in 1:nrow(mean.call.rates)){
    plot(NULL, xlim=c(-2,2),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5, main = paste(dist.bins[i],'-', dist.bins[i+1],'m',sep=' '))
    abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
    abline(v=0, lty=1, col = 'black')
    polygon(c(tseq,rev(tseq)),c(uppers.boot[i,],rev(lowers.boot[i,])),col=alpha(cols[i],.2),border=NA)
    for(j in 1:dim(mean.call.rates.boot)[3]){
      lines(tseq, mean.call.rates.boot[i,,j], col = cols[i], lwd = .2)
    }
    lines(tseq, mean.call.rates[i,], col = cols[i], lwd = 5, type = 'l', )
  }
  
  #Normalized version (to test for peaks and ignore overall up/down shifts across days)
  mean.call.rates.norm <- t(apply(mean.call.rates, 1, function(x){return(x / max(x))}))
  mean.call.rates.boot.norm <- aperm(apply(mean.call.rates.boot, c(1,3), function(x){return(x / max(x))}), c(2,1,3))
  uppers.boot.norm <- apply(mean.call.rates.boot.norm, c(1,2), FUN = function(x){return(quantile(x,0.975,na.rm=T))})
  lowers.boot.norm <- apply(mean.call.rates.boot.norm, c(1,2), FUN = function(x){return(quantile(x,0.025,na.rm=T))})
  
  
  #make the plot
  quartz(height = 6, width = 12)
  par(mfrow=c(1,length(dist.bins)-1))
  par(mar=c(6,6,3,1))
  ymax <- max(mean.call.rates.boot.norm)*1.01
  cols <- viridis(nrow(mean.call.rates.norm))
  for(i in 1:nrow(mean.call.rates.norm)){
    plot(NULL, xlim=c(-2,2),ylim=c(.5,ymax),xlab='Time lag (sec)', ylab = 'Call rate / Max call rate', cex.axis=1.5,cex.lab=1.5, main = paste(dist.bins[i],'-', dist.bins[i+1],'m',sep=' '))
    abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
    abline(v=0, lty=1, col = 'black')
    polygon(c(tseq,rev(tseq)),c(uppers.boot.norm[i,],rev(lowers.boot.norm[i,])),col=alpha(cols[i],.2),border=NA)
    for(j in 1:dim(mean.call.rates.boot.norm)[3]){
      lines(tseq, mean.call.rates.boot.norm[i,,j], col = cols[i], lwd = .2)
    }
    lines(tseq, mean.call.rates.norm[i,], col = cols[i], lwd = 5, type = 'l', )
  }
  
  #---PLOT 3: Call response dynamics across all distances 
  #By age class 
  
  #adult vs adult
  idxs <- which((callresp$caller != callresp$responder) & 
                  (callresp$age_caller %in% adult_classes) & 
                  (callresp$age_responder %in% adult_classes))
  mean.call.rates.adult.adult <- colMeans(callresp.seqs[idxs,],na.rm=T)
  
  #juv vs juv
  idxs <- which((callresp$caller != callresp$responder) & 
                  (callresp$age_caller %in% juv_classes) & 
                  (callresp$age_responder %in% juv_classes))
  mean.call.rates.juv.juv <- colMeans(callresp.seqs[idxs,],na.rm=T)
  
  #adult caller, juv responder
  idxs <- which((callresp$caller != callresp$responder) & 
                  (callresp$age_caller %in% adult_classes) & 
                  (callresp$age_responder %in% juv_classes))
  mean.call.rates.adult.juv <- colMeans(callresp.seqs[idxs,],na.rm=T)

  #juv caller, adult responder
  idxs <- which((callresp$caller != callresp$responder) & 
                  (callresp$age_caller %in% juv_classes) & 
                  (callresp$age_responder %in% adult_classes))
  mean.call.rates.juv.adult <- colMeans(callresp.seqs[idxs,],na.rm=T)
  
  #plot
  quartz(height = 8, width = 12)
  par(mfrow=c(1,2))
  par(mar=c(8,6,3,1))
  ymax_adult <- max(c(mean.call.rates.adult.adult, mean.call.rates.juv.adult))*1.1
  ymax_juv <- max(c(mean.call.rates.juv.juv, mean.call.rates.adult.juv))*1.1
  plot(NULL, xlim=c(-2,2),ylim=c(0,ymax_adult),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5, main = 'All distances')
  abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
  abline(v=0, lty=1, col = 'black')
  lines(tseq, mean.call.rates.adult.adult, col = 'red', lwd = 3, type = 'l', lty = 1)
  lines(tseq, mean.call.rates.juv.adult, col = 'red', lwd = 3, type = 'l', lty = 3)
  legend('bottomleft', col = c('red','red'), lty = c(1, 3), lwd = 3, legend = c('Adult response to Adult', 'Adult response to Juvenile'), bg = 'white')
  
  plot(NULL, xlim=c(-2,2),ylim=c(0,ymax_juv),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5, main = 'All distances')
  abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
  abline(v=0, lty=1, col = 'black')
  lines(tseq, mean.call.rates.juv.juv, col = 'black', lwd = 3, type = 'l', lty = 3)
  lines(tseq, mean.call.rates.adult.juv, col = 'black', lwd = 3, type = 'l', lty = 1)
  legend('bottomleft',col = c('black','black'), lty = c(1, 3), lwd = 3, legend = c('Juvenile response to Adult','Juvenile response to Juvenile'), bg = 'white')
}

#---PLOT 4: Self-reply dynamics (after how long do individuals repeat themselves?) - adult vs juv
if(plot.self.resp.all){
  idxs.adult <- which(callresp$caller == callresp$responder & callresp$age_responder %in% adult_classes)
  idxs.juv <- which(callresp$caller == callresp$responder & callresp$age_responder %in% juv_classes)
  self.reply.rates.adult <- colMeans(callresp.seqs[idxs.adult,], na.rm=T)
  self.reply.rates.juv <- colMeans(callresp.seqs[idxs.juv,], na.rm=T)
  
  #make the plot
  quartz(height = 8, width = 12)
  par(mfrow=c(1,2))
  par(mar=c(8,6,3,1))
  ymax_adult <- max(self.reply.rates.adult)*1.1
  ymax_juv <- max(self.reply.rates.juv)
  
  #adult
  plot(NULL, xlim = c(-5,5), ylim = c(0,max(self.reply.rates.adult)*1.1), lwd = 3, xlab = 'Time (s)', ylab = 'Self-reply rate - Adult', cex.axis = 1.5, cex.lab = 2)
  abline(v = seq(-20,20,1), col = 'gray', lwd = 0.5)
  abline(v=0, lty = 2)
  lines(tseq, self.reply.rates.adult, lwd = 3)
  
  plot(NULL, xlim = c(-5,5), ylim = c(0,max(self.reply.rates.juv)*1.1), lwd = 3, xlab = 'Time (s)', ylab = 'Self-reply rate - Juvenile', cex.axis = 1.5, cex.lab = 2)
  abline(v = seq(-20,20,1), col = 'gray', lwd = 0.5)
  abline(v=0, lty = 2)
  lines(tseq, self.reply.rates.juv, lwd = 3)
  
  
}

#------SUBSET BY CALLER-----
if(plot.by.caller){
  
  #max distance to use
  max.dist <- 5
  
  #get data frame of all callers and their age / dominance classes
  callresp$caller_id_age <- paste(callresp$caller, callresp$age_caller, sep = '_')
  ages <- aggregate(callresp$age_caller, by = list(callresp$caller_id_age), FUN = unique)
  ids <- aggregate(callresp$caller, by = list(callresp$caller_id_age), FUN = unique)
  callers <- data.frame(id_age = ages$Group.1, id = ids$x, age = ages$x)
  
  #get responder call rate as a function of time since caller call
  mean.call.rates <- matrix(NA, nrow = nrow(callers), ncol = length(tseq))
  n.samps <- rep(NA, nrow(callers))
  for(i in 1:nrow(callers)){
    idxs <- which(callresp$caller_id_age==callers$id_age[i] & callresp$distance <= max.dist & callresp$caller != callresp$responder)
    n.samps[i] <- length(idxs)
    mean.call.rates[i,] <- colMeans(callresp.seqs[idxs,])
  }
  
  #make a plot
  quartz(height = 8, width = 12)
  par(mfrow=c(5, 5),mar=c(1,1,1,1))
  for(i in 1:nrow(callers)){
    if(n.samps[i] >= 500){
      plot(tseq, mean.call.rates[i,], type = 'l', lwd = 2, xlim = c(-2,2), xlab = 'Time (sec)', ylab = 'Mean call rate', main = callers$id_age[i])
      abline(v=0)
    }
  }
  
}






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