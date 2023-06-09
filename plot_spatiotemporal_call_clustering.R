#Code to plot the results of spatiotemporal call clustering

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------

savedir <- '~/Dropbox/meerkats/processed_data_serverdownload_2023-01-09/'
testflag <- F
callType <- 'cc'

#----------YOU SHOULD GENERALLY NOT NEED TO MODIFY THESE PARAMETERS--------------
sessions <- c('HM2017','HM2019','L2019')
plot_signif_stars <- F

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
  #first find out when data is significantly greater than or less than null
  signif <- matrix(nrow = nrow(dKdAdt.data), ncol = ncol(dKdAdt.data))
  for(i in 1:nrow(dKdAdt.data)){
    for(j in 1:ncol(dKdAdt.data)){
      data.val <- dKdAdt.data[i,j]
      rand.vals <- dKdAdt.rand[i,j,]
      rand.upper <- quantile(rand.vals, 0.975)
      rand.lower <- quantile(rand.vals, 0.025)
      signif[i,j] <- (data.val > rand.upper) | (data.val < rand.lower) #two tailed
    }
  }
  
  #make the plot
  quartz(height = 8, width = 8)
  par(mfrow=c(1,1), cex.main = 2, cex.lab=2, mar = c(5,5,1,1), cex.axis = 1.5)
  dKdAdt.rand.mean <- apply(dKdAdt.rand, c(1,2), mean)
  ratio <- log(dKdAdt.data / dKdAdt.rand.mean, base = 10)
  ratio[is.infinite(ratio)] <- NA #if either data or randomizations were 0, don't plot (this results in infinite log ratio)
  image.plot(ratio, zlim = c(-max(abs(ratio),na.rm=T),max(abs(ratio),na.rm=T)), col = bluered(256), xlab = 'Distance (m)', ylab = 'Time (sec)', xaxt = 'n', yaxt = 'n', main = session)
  matrix_x_locs <- seq(0,1,length.out = nrow(dKdAdt.data))
  matrix_y_locs <- seq(0,1,length.out = ncol(dKdAdt.data))
  if(plot_signif_stars){
      for(i in 1:length(matrix_x_locs)){
        for(j in 1:length(matrix_y_locs)){
          if(signif[i,j]){
            points(matrix_x_locs[i],matrix_y_locs[j],pch = 8)
          }
        }
      }
  }
  axis(side = 1, at = seq(0,1,length.out = length(dist.windows)), labels = dist.windows)
  axis(side = 2, at = seq(0,1,length.out = length(time.windows)), labels = time.windows)
}
