#Code to plot the results of spatiotemporal call clustering

#-------PARAMETERS------

sessions <- c('HM2017','HM2019','L2019')
savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/'
testflag <- T
callType <- 'cc'

#--------Libraries------
library(gplots)
library(fields)
library(viridis)

#--------------PLOTTING---------

for(sess.idx in 1:length(sessions)){
  
  session <- sessions[sess.idx]
  
  setwd(savedir)
  
  if(testflag){
    load(paste0(savedir,callType, '_clustering_',session,'_test.RData'))
  } else{
    load(paste0(savedir,callType, '_clustering_',session,'.RData'))
  }
  
  #dK/(dAdt) vs distance (panels = time)
  quartz(height = 6, width = 16)
  par(mfrow=c(1,length(time.windows)), mar = c(6,6,2,1))
  for(t in 1:length(time.windows)){
    plot(NULL, xlim = range(dist.windows), ylim = c(0,max(dKdAdt.data)), cex.axis = 1.5, cex.lab = 2, xlab = 'Distance (m)', ylab = 'Clustering of calls (dK / area)', main = paste(session,': dt =', time.windows[t],'sec'), log = 'x')
    abline(v = seq(5,max(dist.windows)+5,5), col = 'gray')
    for(n in 1:n.rands){
      lines(dist.windows, dKdAdt.rand[,t,n], lwd = 0.5, col = 'black')
    }
    lines(dist.windows, dKdAdt.data[,t], lwd=2, col = 'red', pch = 19)
  }
  
  #dKdAdt vs time and distance, log(data / mean(null))
  quartz(height = 8, width = 8)
  par(mfrow=c(1,1), cex.main = 2, cex.lab=2, mar = c(5,5,1,1), cex.axis = 1.5)
  dKdAdt.rand.mean <- apply(dKdAdt.rand, c(1,2), mean)
  ratio <- log(dKdAdt.data / dKdAdt.rand.mean, base = 10)
  image.plot(ratio, zlim = c(-max(abs(ratio),na.rm=T),max(abs(ratio),na.rm=T)), col = bluered(256), xlab = 'Distance (m)', ylab = 'Time (sec)', xaxt = 'n', yaxt = 'n', main = session)
  axis(side = 1, at = seq(0,1,length.out = length(dist.windows)), labels = dist.windows)
  axis(side = 2, at = seq(0,1,length.out = length(time.windows)), labels = time.windows)
}
