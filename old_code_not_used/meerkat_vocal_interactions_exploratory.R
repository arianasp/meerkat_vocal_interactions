#making some exploratory plots of vocal interactions

#libraries
library(RColorBrewer)
library(fitdistrplus)
library(viridis)
library(fields)

#input files
audio.file <- '~/Dropbox/meerkats/meerkats_shared/data/2019_ALL_CALLS_SYNCHED.csv'
gps.file <- '~/Dropbox/meerkats/meerkats_shared/data/HM_COORDINATES_2019_sessions.RData'
ind.file <- '~/Dropbox/meerkats/meerkats_shared/data/HM_INDIVIDUAL_INFO_2019.txt'


#load data
load(gps.file)
calls.all <- read.csv(audio.file, header=T, sep='\t', stringsAsFactors=F)

#PREPROCESS

#parse call types
#add a column for simple call type (parse calls in a sort of hierarchical way - should discuss this in more detail)
calls.all$callSimple <- 'oth'
calls.all$callSimple[which(calls.all$callType %in% c('s','sn','sx','sc','s?'))] <- 's'
calls.all$callSimple[which(grepl('agg',calls.all$callType,ignore.case=T))] <- 'agg'
calls.all$callSimple[which(grepl('so',calls.all$callType,ignore.case=T))] <- 'soc'
calls.all$callSimple[which(grepl('chat',calls.all$callType,ignore.case=T))] <- 'chat'
calls.all$callSimple[which(grepl('mo',calls.all$callType,ignore.case=T))] <- 'mo'
calls.all$callSimple[which(grepl('ld',calls.all$callType,ignore.case=T))] <- 'ld'
calls.all$callSimple[which(grepl('al',calls.all$callType,ignore.case=T))] <- 'al'
calls.all$callSimple[which(grepl('cc',calls.all$callType,ignore.case=T))] <- 'cchyb'
calls.all$callSimple[which(calls.all$callType %in% c('cc','cc*','cc+','ccx','ccc','c'))] <- 'cc'
calls.all$callSimple[which(grepl('sync',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('beep',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('skip',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('start',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('stop',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('digging',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('eating',calls.all$callType,ignore.case=T))] <- NA
calls.all$callSimple[which(grepl('pause',calls.all$callType,ignore.case=T))] <- NA

#convert times to POSIXlt
calls.all$t0 <- as.POSIXlt(calls.all$t0GPS_UTC, tz = 'UTC')
calls.all$tf <- calls.all$t0 + calls.all$duration

#remove double-marked calls (same start and end time and same individual calling)
dups <- which(duplicated(cbind(calls.all$ind,as.numeric(calls.all$t0),as.numeric(calls.all$tf))))
calls.all <- calls.all[-dups,]

#get unique dates
dates <- unique(calls.all$date)

#get a date to work with
date <- dates[5]

#get data from that date
data.date <- calls.all[which(calls.all$date==date),]

#get start and end time
starts <- aggregate(data.date$t0, by = list(data.date$ind), FUN = min, na.rm=T)
ends <- aggregate(data.date$tf, by = list(data.date$ind), FUN = max, na.rm=T)
start.first <- min(starts$x)
start.all <- max(starts$x)
end.all <- min(ends$x)
end.last <- max(ends$x)

#get calls only (no synch calls etc) - NOTE eating and chewing and sync and pause are getting marked as isCall==1
calls <- data.date[which(!data.date$nonFocal & !data.date$unsureFocal & data.date$isCall & !is.na(data.date$callSimple)),]

#colors for simple call cateogries
call.types <- c('mo','agg','soc','al','ld','chat','oth','cchyb','cc','s')
call.colors <- c(brewer.pal(8,'Dark2'),c('black','#0000AA'))


#get unique individuals, give them colors
inds <- unique(calls$ind)
ind.cols <- brewer.pal(length(inds),'Paired')


#CALL TIMESERIES VISUALIZATION PLOT
#make a visualization of calls over time with individuals as rows
quartz(width = 20, height = 6)
par(mar=c(4,12,4,0))
plot(NULL, xlim = c(0, as.numeric(end.last) - as.numeric(start.first)), ylim = c(0,length(inds)), yaxt='n', ylab = NA,xlab='Time (sec)', main = date)
for(i in 1:length(inds)){
  calls.ind <- calls[which(calls$ind == inds[i]),]
  arrows(x0 = as.numeric(calls.ind$t0)-as.numeric(start.first), y0 = i-1, x1 = as.numeric(calls.ind$t0)-as.numeric(start.first), y1 = i, col = call.colors[match(calls.ind$callSimple,call.types)], length = 0, lwd = .5)
}

axis(side = 2, at = seq(.5,length(inds)),labels = inds, las = 1)
abline(h=seq(0,length(inds)))
legend('topright',legend = call.types, col = call.colors, lty = 1, bg = 'white')

#DISTRIBUTION OF INTER-CALL INTERVALS
#look at inter-call intervals for cc's

#get ccs from all dates and call inds
calls.all <- calls.all[order(calls.all$t0),]
ccs <- calls.all[which(!calls.all$nonFocal & !calls.all$unsureFocal & calls.all$isCall & calls.all$callSimple %in% c('cc','cchyb')),]

#create a list of inter-call intervals
icis <- list()
for(i in 1:length(inds)){
  curr <- ccs[which(ccs$ind==inds[i]),]
  icis.curr <- diff(as.numeric(curr$t0))
  icis.curr <- icis.curr[which(icis.curr < 10*60*60)] #remove icis that cut across days
  icis[[i]] <- icis.curr
}

#create plot
quartz(width = 12, height = 8)
#par(mfrow=c(3,4))
plot(NULL,xlim=c(0,60),ylim=c(0,.18),xlab = 'ICI (sec)', ylab = 'Density', cex.lab=1.5)
for(i in 1:length(inds)){
  #lnormfit <- fitdist(data = icis[[i]], distr='lnorm')
  ks <- density(icis[[i]],bw=.5)
  lines(ks$x,ks$y,type='l', col=ind.cols[i] ,lwd=2)
}
legend('topright',legend=inds,col=ind.cols,lty=1,lwd=2)

#are ICIs autocorrelated across time?
quartz(width = 12, height = 8)
par(mfrow=c(3,4))
for(i in 1:length(inds)){
  ici.ind <- icis[[i]]
  acf(ici.ind, main = inds[i])
}

#delay period plots (see Ravignani 2019)

#get data from two individuals (a 'focal' and a 'neighbor')
ind.foc <- inds[1]
ind.nbr <- inds[2]

prev.period.nbr <- c() #holds inter-onset-intervals for neighbor calls
prev.period.foc <- c() #holds inter-onset interval for focal call that was closest in time
delay.foc <- c() #holds delay between neighbor call onset and next call by focal

#loop over neighbor calls to get previous period and next call by focal
for(d in 1:length(dates)){
  date <- dates[d]
  cc.foc <- sort(as.numeric(ccs$t0[which(ccs$ind == ind.foc & ccs$date == date)]))
  cc.nbr <- sort(as.numeric(ccs$t0[which(ccs$ind == ind.nbr & ccs$date == date)]))
  
  #remove time when both were not labeled (NOTE: this is quick and dirty - should fix to start at "START" time instead of first CC)
  min.time <- max(min(cc.foc),min(cc.nbr))
  max.time <- min(max(cc.foc),max(cc.nbr))
  
  #remove data from before min.time or after max.time
  cc.foc <- cc.foc[which(cc.foc >= min.time & cc.foc <= max.time)]
  cc.nbr <- cc.nbr[which(cc.nbr >= min.time & cc.nbr <= max.time)]
  
  #loop over calls and get data for delay period plot
  for(i in 2:length(cc.nbr)){
    t0 <- cc.nbr[i]
    ioi <- t0 - cc.nbr[i-1]
    foc.calls.after.nbr.onset <- which(cc.foc >= t0)
    if(length(foc.calls.after.nbr.onset) > 0){
      t.foc.next <- min(cc.foc[foc.calls.after.nbr.onset])
      
      idx.next.foc.call <- which(cc.foc == t.foc.next)
      
      #add to list
      prev.period.nbr <- c(prev.period.nbr, ioi) 
      delay.foc <- c(delay.foc, t.foc.next - t0)
      
      if(idx.next.foc.call > 1){
        prev.period.foc <- c(prev.period.foc, cc.foc[idx.next.foc.call] - cc.foc[idx.next.foc.call-1])
      } else{
        prev.period.foc <- c(prev.period.foc, NA)
      }
    }
    
  }
  
}

#make the plot
plot(prev.period.nbr, delay.foc, pch = 19, cex = 0.2, xlim=c(0,60),asp=1)
lines(c(0,60),c(0,60))
plot(prev.period.nbr, prev.period.foc, pch = 19, cex = 0.2, xlim=c(0,120),ylim=c(0,120))


#get call rate with exponential decay
decay.rate <- 1.05
date <- dates[1]
ccs.date <- ccs[which(ccs$date==date),]
min.time <- min(as.numeric(ccs.date$t0))
max.time <- max(as.numeric(ccs.date$t0))
tseq <- seq(min.time,max.time,by=1)
call.seqs <- matrix(NA,nrow=length(inds),ncol=length(tseq))
curr.rate <- 0
for(i in 1:length(inds)){
  ccs.ind <- as.numeric(ccs.date$t0[which(ccs.date$ind==inds[i])])
  for(j in 1:(length(tseq)-1)){
    t <- tseq[j]
    tf <- tseq[j+1]
    new.calls <- sum(ccs.ind >= t & ccs.ind < tf)
    curr.rate <- curr.rate / decay.rate + new.calls
    call.seqs[i,j] <- curr.rate
  }
  
}

plot(call.seqs[1,],type='l',xlim=c(5000,6000))
