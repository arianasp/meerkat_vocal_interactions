#Code to plot the results of spatiotemporal call clustering
#From the manuscript "Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups" (working title)
#This code requires that you first run the script get_spatiotemporal_call_clustering.R or use the precomputed output

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------

savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/paper_data_to_submit/precomputed/'
testflag <- F
callType <- 'cc'

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------
sessions <- c('HM2017','HM2019','L2019')
plot_signif_stars <- T

#--------LIBRARIES------
library(gplots)
library(fields)
library(viridis)

#-----------GRAPHICS--------
#make compatible with windows OS
if(.Platform$OS.type=="windows") {
  quartz<-function() windows()
}

#--------------PLOTTING---------

for(sess.idx in 1:length(sessions)){
  
  session <- sessions[sess.idx]
  
  setwd(savedir)
  
  #load data (if testflag == T, load test data, but ignore the results!)
  if(testflag){
    load(paste0(savedir,callType, '_clustering_',session,'_test.RData'))
  } else{
    load(paste0(savedir,callType, '_clustering_',session,'.RData'))
  }
  
  #K vs time and distance (heat map plot) comparing data to null model - heat map shows log(data / mean(null))
  
  #first find out when data is significantly greater than or less than null (alpha = 0.05)
  #signif is a matrix which is T if data significantly different than null and F otherwise
  signif <- matrix(nrow = nrow(K.data), ncol = ncol(K.data))
  for(i in 1:nrow(K.data)){
    for(j in 1:ncol(K.data)){
      data.val <- K.data[i,j]
      rand.vals <- K.rand[i,j,]
      rand.upper <- quantile(rand.vals, 0.975, na.rm=T)
      rand.lower <- quantile(rand.vals, 0.025, na.rm=T)
      signif[i,j] <- (data.val > rand.upper) | (data.val < rand.lower) #two tailed
    }
  }
  
  #make the plot
  quartz(height = 8, width = 8)
  par(mfrow=c(1,1), cex.main = 2, cex.lab=2, mar = c(5,5,1,1), cex.axis = 1.5)
  
  #mean of null model randomizations
  K.rand.mean <- apply(K.rand, c(1,2), mean)
  
  #log ratio between data and null
  ratio <- log(K.data / K.rand.mean, base = 10)
  ratio[is.infinite(ratio)] <- NA #if either data or randomizations were 0, don't plot (this results in infinite log ratio)
  
  #make plot
  image.plot(ratio, zlim = c(-max(abs(ratio),na.rm=T),max(abs(ratio),na.rm=T)), col = bluered(256), xlab = 'Distance (m)', ylab = 'Time (sec)', xaxt = 'n', yaxt = 'n', main = session)
  
  #add significance stars
  matrix_x_locs <- seq(0,1,length.out = nrow(K.data))
  matrix_y_locs <- seq(0,1,length.out = ncol(K.data))
  if(plot_signif_stars){
      for(i in 1:length(matrix_x_locs)){
        for(j in 1:length(matrix_y_locs)){
          if(!is.na(signif[i,j]) & signif[i,j]){
            points(matrix_x_locs[i],matrix_y_locs[j],pch = 8)
          }
        }
      }
  }
  
  #add axes 
  axis(side = 1, at = seq(0,1,length.out = length(dist.windows)), labels = dist.windows)
  axis(side = 2, at = seq(0,1,length.out = length(time.windows)), labels = time.windows)
}
