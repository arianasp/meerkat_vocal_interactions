
setwd('~/Downloads/')

load("tree_data.RData")

#parameters specific to cc or sns
data_list <- list(nodes_cc, nodes_sn)
col_list <- c('red','blue')
titles <- c('Close call sequences','Short note sequences')

#general params
cex <- 7
pch <- 21
scalefac <- 30
switch <- 'gold1'
stay <- 'black'
stop <- 'gray'
linecols <- c(stay, switch, stop, stay, switch, stop, switch, stay, switch)
lty <- rep(1, length(linecols))

#make plots
quartz(width = 12, height = 8)
par(mfrow=c(1,2), mar=c(0,0,0,0))
for(i in 1:length(data_list)){
  #get node data
  nodes <- data_list[[i]]
  col <- col_list[i]
  nodes$p[1] <- sum(nodes$p[3:5])
  nodes$p[2] <- sum(nodes$p[6:9])
  
  #make plots
  plot(NULL, xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(-.1,1.3),ylim=c(-0.1,1.1), bty='n')
  text(.6, 1.1, titles[i], col = col, cex = 2)
  arrows(x0 = c(0,0), y0 = c(.5,.5), x1 = nodes$x[1:2], y1 = nodes$y[1:2], 
         lwd = nodes$p[1:2]*scalefac, length=0,
         col=linecols[1:2])
  arrows(x0 = c(rep(nodes$x[1],3), rep(nodes$x[2],4)), 
         y0 = c(rep(nodes$y[1],3), rep(nodes$y[2],4)),
         x1 = nodes$x[3:9], 
         y1 = nodes$y[3:9], 
         lwd = nodes$p[3:9]*scalefac, length = 0,
         col = linecols[3:9], lty = lty[3:9])
  points(nodes$x, nodes$y, cex = cex, pch = pch, col = 'black', bg=col,lwd=2)
  points(0,.5, cex = cex, pch = pch, col = 'black', bg = col, lwd=2)
  text(nodes$x, nodes$y, labels = nodes$text, col='white',cex=cex.text)
  text(0,.5,'A',col='white',cex=cex.text)
  text(1.15, nodes$y[3:9], paste0(nodes$p[3:9]*100,'%'),cex=cex.text)
  text(nodes$x[1:2], nodes$y[1:2] + c(.07,-.07),paste0(nodes$p[1:2]*100,'%'),cex=cex.text)
}
