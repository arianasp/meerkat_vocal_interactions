#Code to plot the results of spatiotemporal call clustering
#From the manuscript "Mapping vocal interactions in space and time differentiates call-and-response vs. broadcast signalling in meerkat groups" (working title)
#This code requires that you first run the script get_spatiotemporal_call_clustering.R or use the precomputed output

#-------YOU WILL NEED TO MODIFY THESE PARAMETERS TO RUN ON YOUR MACHINE-----------

#savedir <- '~/Dropbox/meerkats/processed_data/vocal_interactions_paper_data_submitted/precomputed/'
savedir <- '~/Desktop/meerkat_data_anon/precomputed/'
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
setwd(savedir)

#Plot by group
for(sess.idx in 1:length(sessions)){
  
  session <- sessions[sess.idx]
  
  #load data (if testflag == T, load test data, but ignore the results!)
  if(testflag){
    load(paste0(savedir,callType, '_clustering_',session,'_test.RData'))
  } else{
    load(paste0(savedir,callType, '_clustering_',session,'_anon.RData'))
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

#PLOT: Aggregate across groups
quartz(height = 8, width = 8)
par(mfrow=c(1,1), cex.main = 2, cex.lab=2, mar = c(8,6,4,1), cex.axis = 1.5)

pairs.agg.data <- matrix(0, nrow = length(dist.windows)-1, ncol = length(time.windows)-1)
pairs.agg.rand <- array(0, dim = c(length(dist.windows)-1, ncol = length(time.windows)-1, n.rands))
tot.pairs.agg.data <- 0
tot.pairs.agg.rand <- rep(0, n.rands)
for(sess.idx in 1:length(sessions)){
  
  session <- sessions[sess.idx]
  
  #load data (if testflag == T, load test data, but ignore the results!)
  if(testflag){
    load(paste0(savedir,callType, '_clustering_',session,'_test.RData'))
  } else{
    load(paste0(savedir,callType, '_clustering_',session,'_anon.RData'))
  }
  
  pairs.agg.data <- pairs.agg.data + K.data*tot.pairs.data
  for(n in 1:n.rands){
    pairs.agg.rand[,,n] <- pairs.agg.rand[,,n] + K.rand[,,n]*tot.pairs.rand[n]
  }
  tot.pairs.agg.data <- tot.pairs.agg.data + tot.pairs.data
  tot.pairs.agg.rand <- tot.pairs.agg.rand + tot.pairs.rand
  
}  

K.data <- pairs.agg.data / tot.pairs.agg.data
K.rand <- array(NA, dim = dim(pairs.agg.rand))
for(n in 1:n.rands){
  K.rand[,,i] <- pairs.agg.rand[,,n] / tot.pairs.agg.rand[n]
}

K.rand.mean <- apply(K.rand, c(1,2), mean, na.rm=T)

colramp <- colorRampPalette(c('darkorchid4','white','darkgoldenrod1'))
cols <- colramp(256)

xaxis <- seq(-1/(length(dist.windows)-1)/2, 1 + 1/(length(dist.windows)-1)/2, length.out = length(dist.windows))
yaxis <- seq(-1/(length(time.windows)-1)/2, 1 + 1/(length(time.windows)-1)/2, length.out = length(time.windows))

mat_to_plot <- log(K.data / K.rand.mean)
mat_to_plot[is.infinite(mat_to_plot)] <- NA
scale <- max(abs(mat_to_plot),na.rm=T)

lab.main <- ''
col.main <- 'black'
if(callType=='cc'){
  lab.main <- 'Close calls'
  col.main <- 'red'
}
if(callType=='sn'){
  lab.main <- 'Short note calls'
  col.main <- 'blue'
}


image.plot(mat_to_plot, col = cols, zlim = c(-scale,scale),xlab = 'Distance (m)', ylab = 'Time (sec)', xaxt = 'n', yaxt = 'n', cex.lab = 1.5, main = lab.main, col.main = col.main)
axis(side = 1, at = xaxis, labels = dist.windows, las = 1, cex = 1.5)
axis(side = 2, at = yaxis, labels = time.windows, las = 1, cex = 1.5)
