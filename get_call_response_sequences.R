#This script is part of an analysis that uses cross-correlation functions to explore call/response dynamics as a function time and space

#The script collects up call-response sequences across all currently available datasets, for subsequent analysis

#It works as follows:
#Make a table of all calls within a given call category (set of calls), identify the one who gave the call as the 'caller'
#For each call, identify all other individuals ('responders') who were present (calls labeled) on that day, note their distance from the caller
#Also identify the distance between the caller and responder at the time when the call was given by the caller
#Save this table as callresp
#Note: This table DOES include self-responses (which can be useful for looking at self-response patterns, or can be filtered out later)

#Then for each of row in this table, look at the call sequence of the 'responder' over time, lined up so that t = 0 = time of the 'caller' call
#Run a kernel over this time series to smooth it with a bandwidth bw = .05 sec (or as specified by parameters bw)
#Save the output to a matrix as callresp.seqs, such that the indices in the table match the indices in the matrix

#In a subsequent script (plot_call_response_sequences), these sequences are further explored and plotted.
#For instance, we can get the mean across all (or a subset) of rows of the matrix for different conditions:
# Different distance bins (aggregate only data from caller-responder pairs within a range of distances from one another)
# Different individuals + distance bins (aggregate only data from caller-responder pairs where the caller was a certain individual, and distance < 3 m since this was discovered to be the relevant range from the distance analysis)
# Investigate self-interactions
# Etc.

#This script requires filling in the parameters specified under PARAMETERS
#If save.output == T, this script saves a file containing:
# callresp: [data frame] of call-response pairs w/ relevant columns
#  date: [str] the date in format yyyymmdd
#  caller: [str] id of the calling individual whose call triggers the sequence (lined up at t = 0)
#  responder: [str] id of the responding individual whose subsequent calls are being monitored
#  t0GPS: [str] time of call onset, synched up to GPS time (using synch calls)
#  tfGPS: [str] time of call offset, synched up to GPS time (using synch calls)
#  t0: [numeric] time in numeric form (seconds since 1970-01-01)
#  tf: [numeric] time in numeric form (seconds since 1970-01-01)
#  start: [numeric] time of the start of the data on that day (all individuals recorded), in numeric form
#  end: [numeric] time of the end of the data on that day (all individuals recorded), in numeric form
#  groupyear: [str] group and year mashed together to indicate which dataset's GPS data to use
#  distance: [numeric] distance between caller and responder at the time of the original call
# callresp.seqs: [matrix] of call-response sequences, 
#   where each row represents the sequence associated with each row in the callresp data frame
#   and each column represents a point in time
# tseq: [vector] giving the time point associated with the columns of callresp.seqs (in seconds)
# params: [named list] of parameters used in the analysis, as described below

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
general.funcs.filename <- 'meerkat_functions.R'

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

#---------------------SETUP-----------------------

#set computer time zone to UTC, to avoid stupid time zone issues
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
    
    #add "spikes" for each responder call
    if(length(dt.foc)>0){
      for(j in 1:length(dt.foc)){
        call.seq <- call.seq + dnorm(tseq, mean = dt.foc[j], sd = bw)
      }
    }
    
    callresp.seqs[i,] <- call.seq
    
  }
}

timestamp()
print('saving output')
if(save.output){
  save(file = paste0(savedir, '/', savename), list = c('callresp','callresp.seqs','tseq', 'params'))
}

print('done')
timestamp()

