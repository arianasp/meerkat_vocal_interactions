#Process audio label files to determine when audio was labeled and when GPS data was recorded (incorporating out of range as well as skipon/skipoff)

#This script generates and saves two matrices for each group-year:
#gpsOn: [n_inds x n_times matrix] boolean, True if GPS data was recorded for that individual at that time
#audioOn: [n_inds x n_times matrix] boolean, True if audio data was labeled for that individual at that time

#The files are saved to the output file labeled:
#groupyear_DATAPRESENCE_all_sessions.RData

#TODO 1: 
#For SOUNDFOC files, currently we assume that GPS data is correct (i.e. not out of range) if not labelled. 
#This means that if there is not a SOUNDFOC file for a given individual, all their GPS data will be included (gpsOn==T).
#Similarly, if oor labeling ends after some point, everything after it will be assumed to be correct GPS data.
#We could consider making the more conservative assumption that if a SOUNDFOC file is not labelled, all the GPS data should be removed.

#TODO 2: 
#This script has not been fully tested!
#CHECK THIS AGAIN to make sure it is performing correctly with no strange edge cases

#------PARAMS------

#directory where gps data is stored for the project
gps.datadir <- '/Volumes/EAS_shared/meerkat/working/processed/movement/'

#use latest audio label table
label.file <- '/Volumes/EAS_shared/meerkat/working/processed/acoustic/resolve_conflicts/foc_calls_resolved.csv'

#groupyears
groupyears <- c('HM2017','HM2019','L2019')

#where to save data ###########HERE#######
savedir <- '/Volumes/EAS_shared/meerkat/working/processed/acoustic'

#------LIBRARIES-----

library(lubridate)

#------MAIN-----

#set computer time zone to UTC, to avoid stupid time zone issues
Sys.setenv(TZ='UTC')

#load labels data
labels.all <- read.csv(label.file,sep=',',header=T)

#get all label files
label.files <- unique(labels.all$csvFileName)

setwd(gps.datadir)

#Load GPS data
for(g in 1:length(groupyears)){
  groupyear <- groupyears[g]
  load(file = paste(groupyear,'_COORDINATES_all_sessions.RData',sep=''))
}

#-----START MAIN PROCESSING-----
for(g in 1:length(groupyears)){

  #get groupyear
  groupyear <- groupyears[g]
  
  print('---------------------------------------')
  print(groupyear)
  print('---------------------------------------')
  
  #extract relevant objects from GPS data
  allX <- eval(as.name(paste(groupyear, 'allX', sep = '_')))
  indInfo <- eval(as.name(paste(groupyear,'indInfo', sep = '_')))
  timeLine <- eval(as.name(paste(groupyears[g],'timeLine', sep = '_')))
  
  #create matrices to hold info on when each sensor is 'on' (i.e. has data)
  audioOn <- matrix(FALSE, nrow = nrow(allX), ncol = ncol(allX))
  gpsOn <- !is.na(allX)
  
  #get all dates in the time line
  dates <- unique(date(timeLine))
  dates_format <- gsub('-','',dates) #put in right format
  
  #get all individuals in the group
  inds <- indInfo$code
  
  #get all label data for that groupyear
  labels.groupyear <- labels.all[which((labels.all$date %in% dates_format) & (labels.all$ind %in% inds)),]
  
  #get all files for that groupyear
  files.groupyear <- unique(labels.groupyear$csvFileName)
  
  #loop over files to fill in where we have data vs not
  for(f in 1:length(files.groupyear)){
    
    #current file
    file <- files.groupyear[f]
    
    #print file
    print(paste('file:', file))
    
    #current labels
    labels.file <- labels.groupyear[which(labels.groupyear$csvFileName == file),]
    
    #get individual index
    ind.code <- labels.file$ind[1]
    ind <- match(ind.code, indInfo$code)
    
    #get starts, stops, startoor, stopoor, skipon, skipoff
    starts.stops.idxs <- which(labels.file$callType %in% c('start','startoor','oorstart','stop','stopoor','oorstop','oor','bir'))
    skipons <- union(grep('skip on', labels.file$callType), grep('skipon', labels.file$callType))
    skipoffs <- union(grep('skip off', labels.file$callType), grep('skipof', labels.file$callType))
    
    #relabel skip on as skipon and skip off as skipoff
    if(length(skipons)>0){
      labels.file$callType[skipons] <- 'skipon'
    }
    if(length(skipoffs)){
      labels.file$callType[skipoffs] <- 'skipoff'
    }

    #gather all indexes
    idxs.to.use <- unique(c(starts.stops.idxs, skipons, skipoffs))
    
    #create table of all starts, stops, oors, birs, skipons, and skipoffs in the file
    starts.stops <- labels.file[idxs.to.use,]
    
    #sort by time
    starts.stops <- starts.stops[order(starts.stops$t0GPS_UTC),]
    
    #add a column for second-level time
    starts.stops$time.sec <- substr(starts.stops$t0GPS_UTC,1,19)
    
    #get time indexes in matrices
    starts.stops$t.idx <- match(starts.stops$time.sec, timeLine)
    
    #start by assuming oor labeling and audio labeling are off
    oorlabel <- F
    audiolabel <- F
    
    #-------COLLAR FILES (audio only, labels don't affect GPS)-------
    if(!grepl('SOUNDFOC',file)){

      audio <- F
      
      #loop over labels in file
      for(i in 1:nrow(starts.stops)){
        
        #get label and time idx
        label <- starts.stops$callType[i]
        t.idx <- starts.stops$t.idx[i]
        
        if(label %in% c('start','skipoff','skipof')){
          if(!audio){
            tprev <- t.idx
          } #else{
            #warning(paste('collar: encountered start label when audio already on: file',file, ', label', label, 'at',t.idx))
          #}
          audio <- T
        }
        
        if(label %in% c('stop','end','skipon')){
          if(audio){
            if(!is.na(tprev) && !is.na(t.idx) && tprev < (t.idx-1)){
              audioOn[ind, tprev:(t.idx-1)] <- T
            }
          } #else{
            #warning(paste('collar: encountered stop label when audio already off:',label, 'at',t.idx))
          #}
          audio <- F
        }
      }
    }
    
    #--------FOCAL FOLLOW FILES (labels affect both GPS and audio)--------
    if(grepl('SOUNDFOC',file)){
      
      #loop over labels and turn audio and gps on and off accordingly
      for(i in 1:nrow(starts.stops)){
        
        label <- starts.stops$callType[i]
        t.idx <- starts.stops$t.idx[i]
        
        #---start marker means both oor and audio labeling are happening---
        if(label == 'start'){
          oorlabel <- T
          audiolabel <- T
          audio <- T
          gps <- T
          tprev <- t.idx
        }
        
        #---startoor / oorstart turns only gps labeling on, audio still off---
        if(label %in% c('startoor', 'oorstart')){
          oorlabel <- T
          gps <- T
          tprev <- t.idx
        }
        
        #---stop turns audio off, and gps off if there is no later "stopoor" label in the file---
        if(label == 'stop'){
          
          #turn audio labeling off
          audiolabel <- F
          
          #gps: if there are no other stopoor labels later, then quit oor labeling and end gps, otherwise leave it
          if(sum(grepl('stopoor', labels.file$callType[i:nrow(labels.file)])) + sum(grepl('oorstop',labels.file$callType[i:nrow(labels.file)])) == 0){
            oorlabel <- F
            gps <- F
          }
          
          #audio: if audio was previously on, fill in audio until the stop
          if(audio){
            if(!is.na(tprev) && !is.na(t.idx) && tprev < (t.idx-1)){
              audioOn[ind, tprev:(t.idx-1)] <- T
            }
          }
          
          #update t.idx
          tprev <- t.idx
          
          #turn audio off
          audio <- F
          
        }
        
        #---stopoor turns gps off and oorlabel off
        if(label %in% c('stopoor', 'oorstop')){
          oorlabel <- F
          gps <- F
        }
        
        #---out of range---
        if(label == 'oor'){
          
          #if oor labeling is not on, throw warning, turn oorlabel on, then assume oor labeling from that point on
          if(!oorlabel){
            #warning(paste('soundfoc: found oor label when oor labeling was not on, file:', file))
            oorlabel <- T
          }
          
          #gps turns off
          gps <- F
          
          #if audio is also being labeled, deal with audio
          if(audiolabel){
            #if the audio was on, turn it off and record the last segment as on
            if(audio){
              if(!is.na(tprev) && !is.na(t.idx) && tprev < (t.idx-1)){
                audioOn[ind, tprev:(t.idx-1)] <- T
              }
              audio <- F
            }
          }
          
          #set tprev to t.idx
          tprev <- t.idx
        }
        
        #---back in range---
        if(label == 'bir'){
          
          #if oor labeling is not on, throw warning, then turn oorlabel on from that point on
          if(!oorlabel){
            #warning(paste('soundfoc: found bir label when oor labeling was not on, file:', file))
            oorlabel <- T
          } 
          
          #if gps was off, remove previous interval and turn back on
          if(!gps){
            if(!is.na(tprev) && !is.na(t.idx) && tprev < (t.idx-1)){
              gpsOn[ind, tprev:(t.idx-1)] <- F
            }
            gps <- T
          }
          
          #if audio labeling was on, also turn audio on
          if(audiolabel){
            audio <- T
          }
          
          #set tprev to t.idx
          tprev <- t.idx
          
        }
        
        #---skipon---
        #affects audio but not gps
        if(label == 'skipon'){
          
          #if skipon is encountered when audio labeling was supposedly off, throw warning
          if(!audiolabel){
            warning(paste('soundfoc: skipon label encountered when audio labeling was not on in general, file:', file, 'label:', label, 't.idx',t.idx))
          }
          #if(!audio){
          #  warning(paste('soundfoc: skipon label encountered when audio labeling was currently off, file:', file, 'label:', label, 't.idx', t.idx))
          #}
          
          #record audio segment and turn audio off
          if(audio){
            if(!is.na(tprev) &&  !is.na(t.idx) && tprev < (t.idx-1)){
              audioOn[ind, tprev:(t.idx-1)] <- T
            }
          }
          audio <- F
          tprev <- t.idx
  
        }
        
        #---skipoff---
        #affects audio only, not gps
        if(label == 'skipoff'){
          
          #if skipoff is encountered when audio labeling was supposedly off, throw warning
          if(!audiolabel){
            #warning(paste('soundfoc: skipoff label encountered when audio labeling was not on in general, file:', file))
            audiolabel <- T
          }
          
          #if audio was already on, throw warning
          #if(audio){
            #warning(paste('skipoff label encountered when audio labeling was already on, file:', file))
          #}
          
          #if audio was previously off, then turn audio on
          audio <- T
        }
      }
    }
  }
  
  #assign to named variables
  assign(x = paste(groupyear,'gpsOn', sep = "_"), value = gpsOn)
  assign(x = paste(groupyear,'audioOn', sep = "_"), value = audioOn)
  
}

#save
for(g in 1:length(groupyears)){
  groupyear <- groupyears[g]
  save(file = paste(savedir,'/', groupyear, '_DATAPRESENCE_all_sessions.RData',sep=''), list = c(paste(groupyear,'audioOn',sep='_'), paste(groupyear,'gpsOn',sep='_')))
}
#
  