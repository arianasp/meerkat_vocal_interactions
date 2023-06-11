#This script is part of an analysis that uses cross-correlation functions to explore call/response dynamics as a function time and space

#The script collects up call-response sequences across all currently available datasets, for subsequent analysis

#It works as follows:
#Make a table of all calls within a given call category (set of calls), identify the one who gave the call as the 'caller'
#For each call, identify all other individuals ('responders') who were present (calls labeled) on that day
#Also identify the distance between the caller and responder at the time when the call was given by the caller
#Save this table as callresp
#Note: This table DOES include self-responses (which can be useful for looking at self-response patterns, or can be filtered out later)

#Then for each of row in this table, look at the call sequence of the 'responder' over time, lined up so that t = 0 = time of the 'caller' call
#Run a kernel over this time series to smooth it with a bandwidth bw = .1 sec (or as specified by parameters bw)
#Save the output to a matrix as callresp.seqs, such that the indices in the table match the indices in the matrix

#In a subsequent script (plot_call_response_sequences.R), these sequences are further explored and plotted.
#For instance, we can get the mean across all (or a subset) of rows of the matrix for different conditions:
# Different distance bins (aggregate only data from caller-responder pairs within a range of distances from one another)
# Different age/sex/dominance classes (aggregate only data from caller-responder pairs where the caller was a certain status, and the responder was a certain status)
# Investigate self-interactions
# Etc.

#This script requires filling in the parameters specified under PARAMETERS 
#Generally, you will only need to modify the parameters specifying the paths to your input and output files, 
#and the caller.calltype and responder.calltype (should be either 'cc' for close calls or 'sn' for short notes, for both variables)
#This is specified in the section below titled YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE

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

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------
#flag that allows you to test the code on a smaller subset of dates (to save time)
#set to T for testing mode, set to F to run code on all data (takes several hours)
#Note that if you run this code in test mode, the output will be very noisy because only a few calls are used, therefore the results are NOT expected to match those in the paper
testflag <- F

#directory where gps data is stored for the project
gps.datadir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 

#direcotry where audio labeling data is stored for the project
audio.datadir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/' 

#directory where code is stored for this project
codedir <- '~/Dropbox/code_ari/meerkat_vocal_interactions'

#directory of where to save results (for later plotting)
savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/'

#list of call types to include in the set of calls by the initial caller (which determines the 0 point of the correlogram)
#Options are either 'cc' (for close calls) or 'sn' (for short note calls)
caller.calltype <- 'cc' 

#list of call types of include in the set of calls by the responder (determines the curve in the correlogram)
#Options are either 'cc' (for close calls) or 'sn' (for short note calls)
responder.calltype <- 'cc'

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------
#whether to save output
save.output <- T

#bandwidth of smoothing kernel (default 0.1)
bw <- .1

#maximum time lag to consider (time since another individual called)
max.lag <- 3

#time step to use for the sequence of times
step <- .02

#what to trigger on (call begin or call end)
trigger.on <- 'begin'

#list of group years
groupyears <- c('HM2017', 'HM2019', 'L2019')

#file names for the audio file and the gps files for each session
gps.files <- paste(groupyears, 'COORDINATES_all_sessions.RData', sep = '_')
audio.file <- 'all_calls_sync_resolved_with_oor_2022-12-04.csv' 

# ----- SETUP -------
#store parameters in a named list
params <- list()
params$gps.datadir <- gps.datadir
params$codedir <- codedir
params$bw <- bw
params$max.lag <- max.lag
params$step <- step
params$caller.calltype <- caller.calltype
params$responder.calltype <- responder.calltype
params$trigger.on <- trigger.on
params$groupyears <- groupyears
params$audio.file <- audio.file
params$gps.files <- gps.files
params$testflag <- testflag

#set computer time zone to UTC, to avoid stupid time zone issues
Sys.setenv(TZ='UTC')

#if saving output, generate the filename to save as
#if running in test mode, use the _test suffix
if(save.output){
  if(testflag){
    savename <- paste0('callresp_', caller.calltype, '_', responder.calltype, '_bw', bw, '_test.RData')
  } else{
    savename <- paste0('callresp_', caller.calltype, '_', responder.calltype, '_bw', bw, '.RData')
  }
}

#-------LIBRARY-------
library(lubridate)

#------------------LOAD DATA-----------------------

timestamp()
print('loading data')

#load audio data and datapresence data (which is in the audio data folder, because it relies on the audio labels)
setwd(audio.datadir)

#the audio data is stored in a data frame called calls.all
calls.all <- read.csv(audio.file, header=T, sep=',', stringsAsFactors=F)

#data presence data for each session
#gpsOn and audioOn are matrices of size n_inds x n_times
#gpsOn[i,t] is T if the gps was recording for that individual at that time index (otherwise F), similar for audioOn
for(g in 1:length(groupyears)){
  load(paste(groupyears[g],'DATAPRESENCE_all_sessions.RData',sep='_'))
}

#load gps data for all sessions
setwd(gps.datadir)
for(i in 1:length(gps.files)){
  load(gps.files[i])
}

#the gps data loaded consists of 
# -- allX and allY matrices (allX[i,t] = easting position of individual i at time index t), similar for allY
# -- timeLine: a vector of timestamps, associated with the column indexes of the allX and allY data
# --indInfo: a data frame with info about the individual names (code) and age/sex/dominance

#--------------------PREPROCESS-------------------------------

timestamp()
print('preprocessing')

#make sure trigger.on is specified as either triggering on beginning or end of call, otherwise error
if(!(trigger.on %in% c('begin','end'))){
  stop('trigger.on parameter should be either begin or end')
}

#convert times to POSIXlt
calls.all$t0 <- as.POSIXlt(calls.all$t0GPS_UTC, tz = 'UTC')
calls.all$tf <- as.POSIXlt(calls.all$tendGPS_UTC, tz = 'UTC')

#remove double-marked calls (same start and end time and same individual calling)
dups <- which(duplicated(cbind(calls.all$ind,as.numeric(calls.all$t0),as.numeric(calls.all$tf))))
if(length(dups)>0){
  calls.all <- calls.all[-dups,]
}

#get unique dates
dates <- unique(calls.all$date)

#indicate which calls are in the caller call type list, and which are in repsonder call type list
calls.all$isCallerCallType <- F
calls.all$isResponderCallType <- F
calls.all$isCallerCallType[grepl(caller.calltype, calls.all$callType)] <- T
calls.all$isResponderCallType[grepl(responder.calltype, calls.all$callType)] <- T
print('For the caller, using all calls fitting the regular expression:')
print(caller.calltype)
print('For the responder, using all calls fitting the regular expression:')
print(responder.calltype)

#---------------------------MAIN-------------------------------
#time sequence bins - a series of time lags relative to the caller's call (at t=0), ranging from - max.lag to max.lag in steps of step
tseq <- seq(-max.lag, max.lag, step)

#get together a table of all calls we are using in the analysis t0 and tf times from all individuals, as well as date and the start and end time of the labeled sequence
calls.use <- data.frame()

#loop over all collaring sessions (groupyears) to build calls.use table
for(g in 1:length(groupyears)){
  
  #get groupyear
  groupyear <- groupyears[g]
  
  #get dates associated with that groupyear
  timeLine <- eval(as.name(paste(groupyears[g],'timeLine', sep = '_')))
  indInfo <- eval(as.name(paste(groupyear,'indInfo', sep = '_')))
  dates <- unique(date(timeLine))
  
  dates_format <- gsub('-','',dates) #put in right format to match calls data frame
  
  #get all individuals in the group
  inds <- indInfo$code
  
  #for each date
  for(d in 1:length(dates)){
  
    date <- dates_format[d]
    
    #get calls on that date
    calls.date <-  calls.all[which((calls.all$date == date) & (calls.all$ind %in% inds)),]
    
    #make a simple table that just contains the cc's from each individual at that date and their t0 and tf times (as numeric)
    calls.use.date <- calls.date[which((calls.date$isCallerCallType | calls.date$isResponderCallType) & as.character(calls.date$pred_focalType) == 'F'),]

    #skip if there were no calls on that date, otherwise...
    if(nrow(calls.use.date)>0){
      
      #add columns for the t.idx (for linking with audioON table) and groupyear
      calls.use.date$t.sec <- substr(calls.use.date$t0GPS_UTC, 1, 19)
      calls.use.date$t.idx <- match(calls.use.date$t.sec, timeLine)
      calls.use.date$groupyear <- groupyear
      
      #append to the larger table
      calls.use.date <- calls.use.date[,c('date','ind','callType','isCallerCallType','isResponderCallType','t0GPS_UTC','tendGPS_UTC','t0','tf','t.idx','groupyear')]
      calls.use.date$t0 <- as.numeric(calls.use.date$t0)
      calls.use.date$tf <- as.numeric(calls.use.date$tf)
      calls.use <- rbind(calls.use, calls.use.date)
    }
  }
}

#rename column names of the calls.use table for clarity
colnames(calls.use) <- c('date','caller','callType','isCallerCallType','isResponderCallType','t0GPS','tfGPS','t0','tf','t.idx','groupyear')

#if using in test mode, subsample calls randomly to select 1000 calls
if(testflag){
  calls.use <- calls.use[sample(1:nrow(calls.use), 1000),]
}

#create an expanded data frame (callresp) that contains rows for each other 'responder' individual (not the caller)
#so for each call by a given individual, there will be a row for every other individual in the group
timestamp()
print('generating call-response table')
callresp <- data.frame()
for(g in 1:length(groupyears)){
  
  groupyear <- groupyears[g]
  
  #get data and rename to be called audioOn as opposed to {groupyear}_audioON (and same for other variables)
  audioOn <- eval(as.name(paste(groupyear,'audioOn', sep = '_')))
  indInfo <- eval(as.name(paste(groupyear,'indInfo', sep = '_')))
  
  #get only rows that are associated with that groupyear
  rows <- which(calls.use$groupyear == groupyear)
  
  print('processing callresp for groupyear:')
  print(groupyear)
  
  #for each row in the call response table in the given groupyear
  for(i in rows){
    
    row <- calls.use[i,]
    t.idx <- row$t.idx
    
    #if time window falls outside of overall recording period bounds, skip
    n.times <- ncol(audioOn)
    if(is.na(t.idx) || (t.idx - max.lag) < 1 || (t.idx + max.lag) > n.times){
      next
    }
    
    #if the caller is of the type we are using
    if(row$isCallerCallType == T){
      
      #find which individuals were possible responders (they were not missing audio data around that call)
      missingAudio <- rowSums(!audioOn[,(t.idx - max.lag):(t.idx + max.lag)])
      inds <- indInfo$code[which(missingAudio == 0)]
      
      #create rows for each possible (non-missing) responder at that time
      new.rows <- row[rep(1,length(inds)),]
      new.rows$responder <- inds
      
      #add to the bottom of the callresp data frame we are building up
      callresp <- rbind(callresp, new.rows)
    }
    
  }
}
#now we have a callresp data frame with callers and possible repsonders for each call of the correct type given by the correct status class
#now we need to add some other info to this data frame

#get distance between caller and responder at time of the call (t0)
#loop over group years
timestamp()
print('computing distances between individuals across all recording periods')

callresp$distance <- callresp$caller.idx <- callresp$responder.idx <- NA
for(g in 1:length(groupyears)){
  
  #get indices associated with that group year, and associated timeline
  idxs <- which(callresp$groupyear == groupyears[g])

  #get the gps data and indInfo data associated with that groupyear
  allX <- eval(as.name(paste(groupyears[g], 'allX', sep = '_')))
  allY <- eval(as.name(paste(groupyears[g], 'allY', sep = '_')))
  indInfo <- eval(as.name(paste(groupyears[g],'indInfo', sep = '_')))
  gpsOn <- eval(as.name(paste(groupyears[g],'gpsOn', sep = '_')))
  
  #replace the times when the gps data was not valid (e.g. when focal follower was out of range) with NAs
  allX[!gpsOn] <- NA
  allY[!gpsOn] <- NA
  
  #get indices of caller and responder (what rows to access in the spatial data, which match the rows in the indInfo table)
  callresp$caller.idx[idxs] <- match(callresp$caller[idxs], indInfo$code)
  callresp$responder.idx[idxs] <- match(callresp$responder[idxs], indInfo$code)
  
  #compute distance between caller and responder at the time of the call
  xc <- allX[cbind(callresp$caller.idx[idxs], callresp$t.idx[idxs])]
  yc <- allY[cbind(callresp$caller.idx[idxs], callresp$t.idx[idxs])]
  xr <- allX[cbind(callresp$responder.idx[idxs], callresp$t.idx[idxs])]
  yr <- allY[cbind(callresp$responder.idx[idxs], callresp$t.idx[idxs])]
  dists <- sqrt((xc - xr)^2 + (yc - yr)^2)
  
  #store distances in callresp table
  callresp$distance[idxs] <- dists
  
}

#now that we have the matrix of all callers and responders, call times, and distances between them, we need to generate vectors of the vocal sequences of the responders (callresp.seqs)

#create a matrix to hold call/response sequences
#each row will represent a sequence of a particular responder associated with the call of a caller
#the rows match up with the callresp data frame
callresp.seqs <- matrix(data=NA,nrow=nrow(callresp),ncol = length(tseq))

#loop through and get individual call patterns starting at tf
timestamp()
print('gathering vocal sequences across all recording periods')
for(i in 1:nrow(callresp.seqs)){
  
  #get the time point that we will consider to be 0 (the time of caller's call, either beginning or end of it)
  if(trigger.on == 'begin'){
    zero.time <- callresp$t0[i]
  } else{
    zero.time <- callresp$tf[i]
  }
  
  #get caller and responder id
  responder <- callresp$responder[i]
  caller <- callresp$caller[i]
  
  #get call times of the responder (only of the correct call type, and removing the current call)
  current.call.idx <- which(calls.use$caller == caller & calls.use$t0 == callresp$t0[i] & calls.use$tf == callresp$tf[i])
  
  #some error checking - if can't find the call throw an error
  if(length(current.call.idx)==0){
    stop('could not find current call in the table')
  } 
  #this should also not happen because we removed duplicates before, but just as a check
  if(length(current.call.idx) > 1){
    warning('found two calls with same label, start, and stop time')
  }
  
  #remove the current call from the table of calls we are using
  calls.use.removeocurrent <- calls.use[-current.call.idx,]
  
  #get call times of the current responder (foc), if they are of the correct type
  #this line might be slightly confusing, but basically we are now looking for the rows where the calls were given by the specified responder (caller == responder)
  #in some cases, the caller and responder might be the same individual - we keep these in the table to be able to look at self-replies later if we want to
  call.times.foc <- calls.use.removeocurrent$t0[which(calls.use.removeocurrent$isResponderCallType & calls.use.removeocurrent$caller == responder)]
  
  #get time differences between the responder call times and the zero time
  dt.foc <- call.times.foc - zero.time
  dt.foc <- dt.foc[which(abs(dt.foc) < max.lag)] #get rid of way early or way late calls to reduce compute time
  
  #generate a vector to hold the call sequence of the responder
  #it will be 0 where there are no calls, with little "spikes" where calls are (we use gaussian kernels with width bw)
  call.seq <- rep(0, length(tseq))
  
  #add "spikes" for each responder call
  if(length(dt.foc)>0){
    for(j in 1:length(dt.foc)){
      call.seq <- call.seq + dnorm(tseq, mean = dt.foc[j], sd = bw)
    }
  }
  
  #add the current sequence to the big matrix of all sequences
  callresp.seqs[i,] <- call.seq
}

timestamp()
print('saving output')
if(save.output){
  save(file = paste0(savedir, '/', savename), list = c('callresp','callresp.seqs','tseq', 'params'))
}

print('done')
timestamp()

