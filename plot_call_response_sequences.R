#This script plots the call-response curves for the paper "Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups"
#BEFORE you run this script, please run the script get_call_response_sequences.R
#This prior script will output a results file calls callresp_xx_xx_bw0.1.RData which contains data required for the plotting
#(xx will either be cc or sn depending on which call type you have selected)
#You will need to specify the path to this results file below.

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------

#file to use for the plots (outputted from get_call_response_sequences)
filename <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/precomputed/callresp_sn_sn_bw0.1.RData'

#data directory where ind info (and gps data) is stored 
datadir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/'

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------
#Distance bins to use for plotting call-response sequences as a function of space
dist.bins <- c(0,2,5,10,50)

n.boots <- 100 #number of bootstraps to do for error bars
gy <- NULL #which group year to use (if NULL, use all groups together)

#whether to plot some extra plots that are not in the paper (default false)
extra.plots <- F

#caller age classes
adult_classes <- c('DominantF','DominantM','Yearling','Sub-Adult','Adult')
juv_classes <- c('Juvenile')

#------------LIBRARIES-------------
library(viridis)

#-----------GRAPHICS--------
#make compatible with windows OS
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

#--------------------HELPER FUNCS---------------
#Calculate mean call rates from data frame of call response data and time sequences
#INPUTS:
# callresp: data frame containing info on call responses (from get_call_response_sequences script)
# callresp.seqs: matrix containing the actual time series data (from get_call_response_sequences script)
# gy: groupyear, currently 'HM2017','HM2019', and 'L2019' are supported. If NULL, use all groups together
# tseq: vector sequence of times (from get_call_response_sequences script)
# dist.bins: vector of distance bins (default to 0,2,5,10,50)
# adult_classes: age classes to use (default to c('DominantF','DominantM','Yearling','Sub-Adult','Adult'))
#OUTPUTS:
# mean.call.rates: matrix of the mean call rate for each distance bin (rows) and time bin (columns)
get_mean_call_rates <- function(callresp, callresp.seqs, gy, tseq, 
                                dist.bins = c(0,2,5,10,50), 
                                adult_classes = c('DominantF','DominantM','Yearling','Sub-Adult','Adult')){
  
  #loop over distance bins to get mean call rates
  mean.call.rates <- matrix(NA,nrow=length(dist.bins)-1, ncol = length(tseq))
  
  for(i in 2:length(dist.bins)){
    
    #get indexes to the rows in the callresp data frame that are within that distance bin, 
    #where the callers and responder are different individuals, and where the caller and responder are both adults
    #if a group year is specified, also include rows only from that groupyear
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
    
    #store mean call rates
    mean.call.rates[i-1,] <- colMeans(callresp.seqs[idxs,],na.rm=T)
  }
  
  #return mean call rates
  return(mean.call.rates)
}

#--------------------LOAD AND PREPROCSS DATA------------------

#load the data computed from the script get_call_response_sequences.R
load(filename)

#get individual info from all groupyears into a single table and add groupyear_ind column
setwd(datadir)
groupyears <- unique(callresp$groupyear)
ind_info <- data.frame()
for(g in 1:length(groupyears)){
  groupyear <- groupyears[g]
  load(paste0(groupyear,'_COORDINATES_all_sessions.RData'))
  ind_info_curr <- eval(as.name(paste(groupyear,'indInfo', sep = '_')))
  ind_info_curr$groupyear_ind <- paste(groupyear, ind_info_curr$code,sep= '_')
  ind_info <- rbind(ind_info, ind_info_curr)
  
  #remove gps data (not needed for further analyses)
  rm(list = paste0(groupyear,c('_allX','_allY','_indInfo','_movSummary','_timeLine')))
}


#add columns to specify age class of all individuals (based on the big ind_info table generated above)
callresp$groupyear_caller <- paste(callresp$groupyear,callresp$caller, sep='_')
callresp$groupyear_responder <- paste(callresp$groupyear, callresp$responder, sep = '_')
callresp$age_caller <- ind_info$status[match(callresp$groupyear_caller, ind_info$groupyear_ind)]
callresp$age_responder <- ind_info$status[match(callresp$groupyear_responder, ind_info$groupyear_ind)]

#---------------------PLOTTING-------------------

#CALL RESPONSE PLOTS
  
#---PLOT 1: Call response dynamics on average, between two adult individuals
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
  
  #bootstrap for error bars - sample days randomly (keep days together because non-independent)
  groupyeardates.boot <- sample(groupyeardates, replace=T)
  
  #get the associated indexes (row in the data frame) associated with the bootstrapped dates
  groupyeardates.boot.idxs <- c()
  for(bb in 1:length(groupyeardates.boot)){
    groupyeardates.boot.idxs <- c(groupyeardates.boot.idxs, which(callresp$groupyeardate == groupyeardates.boot[bb]))
  }
  #make a temporary callresp table and callresp.seqs matrix with the bootstrapped samples
  callresptmp <- callresp[groupyeardates.boot.idxs,]
  callrespseqstmp <- callresp.seqs[groupyeardates.boot.idxs,]
  
  #compute the mean call rates for the bootstrapped sample
  mean.call.rates.boot[,,b] <- get_mean_call_rates(callresptmp, callrespseqstmp, gy, tseq, dist.bins, adult_classes)
}

#get upper and lowers (95% interval) for the bootstrapped samples - actually not really needed as we no longer plot this
uppers.boot <- apply(mean.call.rates.boot, c(1,2), FUN = function(x){return(quantile(x,0.975,na.rm=T))})
lowers.boot <- apply(mean.call.rates.boot, c(1,2), FUN = function(x){return(quantile(x,0.025,na.rm=T))})

#make the plot
quartz(height = 6, width = 12)
par(mfrow=c(1,length(dist.bins)-1))
par(mar=c(6,6,3,1))
ymax <- max(mean.call.rates)*1.3
cols <- viridis(nrow(mean.call.rates))
for(i in 1:nrow(mean.call.rates)){
  
  #set up plot for a given distance range
  plot(NULL, xlim=c(-2,2),ylim=c(0,ymax),xlab='Time lag (sec)', ylab = 'Call rate', cex.axis=1.5,cex.lab=1.5, main = paste(dist.bins[i],'-', dist.bins[i+1],'m',sep=' '))
  abline(v = seq(-3,3,.1), col = 'gray', lwd = 0.5)
  abline(v=0, lty=1, col = 'black')
  
  #plot bootstrapped curves
  for(j in 1:dim(mean.call.rates.boot)[3]){
    lines(tseq, mean.call.rates.boot[i,,j], col = cols[i], lwd = .2)
  }

  #plot full data curves
  lines(tseq, mean.call.rates[i,], col = cols[i], lwd = 5, type = 'l', )
}

#---PLOT 2: Call response dynamics across all distances by age class 

#get indexes associated with different caller/responder age class combos, and compute mean call rates

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
  
# ----- EXTRA PLOTS -----  
#some extra plots that are not included in the paper 
  
if(extra.plots){
#---Spiking neurons plot
  i <- 2
  n.neurons <- 100
  xmax <- 3
  isna <- is.na(rowSums(callresp.seqs))
  idxs <- which((callresp$distance >= dist.bins[i-1]) & (callresp$distance < dist.bins[i]) & !isna)
  idxs.sub <- sample(idxs, n.neurons)
  
  quartz(height = 8, width = 4)
  plot(NULL, xlim = c(-xmax, xmax), ylim = c(0, n.neurons), xlab = "Time (sec)", ylab = 'Calls')
  for(j in 1:n.neurons){
    lines(tseq, callresp.seqs[idxs.sub[j],] + j )
  }
  abline(v = 0, col = 'blue')

  #---Self-reply dynamics (after how long do individuals repeat themselves?) - adult vs juv
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

  #------SUBSET BY CALLER-----
  
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




