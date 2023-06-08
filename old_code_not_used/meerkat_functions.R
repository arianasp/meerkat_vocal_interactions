#Various functions that are useful for analysis of meerkat data

library(circular)
library(dismo)
library(lubridate)
library(gplots)
library(viridis)
library(raster)
library(rgdal)


#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	EastingsCol1: whether eastings should be given in first column of output (default) or not
#Outputs:
#	EastNorths or NorthEasts: [N x 2 matrix] of Eastings and Northings - eastings are first column by default
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='34',southern_hemisphere=TRUE){
	latlons <- data.frame(X=LonsLats[,2],Y=LonsLats[,1])
	non.na.idxs <- which(!is.na(latlons$X) & !is.na(latlons$Y))
	len <- nrow(latlons)
	non.na.latlons <- latlons[non.na.idxs,]
	coordinates(non.na.latlons) <- ~Y + X
	proj4string(non.na.latlons) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')
	if(southern_hemisphere){
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
	} else{
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
	}
	utm <- spTransform(non.na.latlons,CRS(projection.string))
	EastNorths <- matrix(NA,nrow=len,ncol=2)
	EastNorths[non.na.idxs,] <- utm@coords
	if(!EastingsCol1){
		NorthEasts <- EastNorths[,c(2,1)]
		return(NorthEasts)
	} else{
		return(EastNorths)
	}
}

#Converts a matrix of eastings and northings (eastings first column, northings second column) to UTM
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	LonsCol1: whether lons should be given in first column of output (default) or not
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column by default 
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '34',southern_hemisphere=TRUE){
	utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
	non.na.idxs <- which(!is.na(utms$X) & !is.na(utms$Y))
	len <- nrow(utms)
	if(southern_hemisphere){
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
	} else{
		projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
	}
	non.na.utms <- SpatialPoints(utms[non.na.idxs,],proj4string=CRS(projection.string))
	lonlat <- spTransform(non.na.utms,CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
	LonLats <- matrix(NA,nrow=len,ncol=2)
	LonLats[non.na.idxs,] <- lonlat@coords
	if(!LonsCol1){
		LatLons <- LonLats[,c(2,1)]
		return(LatLons)
	} else{
		return(LonLats)
	}	
}

#VISUALIZATIONS

#Function to plot all data from a given burrow (morning, evening, or both), and can also specify a group
#Inputs:
#	burrow: the burrow ID number
#	group: the group ID number (or NULL, if plotting data from all groups)
#	time.of.day: string specifying 'morn','eve', or 'both' (which time(s) of day to plot- defaults to both)
#	burrow.use.over.time: data frame of burrow data over time
#	gps: data frame of GPS data
#	morning.cols, evening.cols: vectors of colors to use for morning and evening data
#	burrow.col: what color to use for the burrow
#	morning.evening.thresh: dividing hour between morning and evening (defaults to 13, meaning 13:00)
#Outpus:
#	A plot of trajectories to/from the burrow
plot.burrow.trajectories <-function(burrow,group=NULL,time.of.day='both',burrow.use.over.time=burrow.use.over.time,gps=gps,morning.cols=adjustcolor(rich.colors(24*60),alpha=0.5),evening.cols=adjustcolor(rich.colors(24*60),alpha=0.5),burrow.col='blue',morning.evening.thresh=13){
	
	#whether to plot morning or evening data or both (default)
	morn <- eve <- T
	if(time.of.day %in% c('morning','morn','Morn','Morning')){
		eve <- F
	}
	if(time.of.day %in% c('evening','eve','Eve','Evening')){
		morn <- F
	}
	
	#get data for the specified burrow (and group, if needed)
	if(!is.null(group)){
		curr.burrow.dat <- burrow.use.over.time[which(burrow.use.over.time$burrow==burrow & burrow.use.over.time$group==group),] 
	} else{
		curr.burrow.dat <- burrow.use.over.time[which(burrow.use.over.time$burrow==burrow),]
	}

	#get dates that burrow was used (evenings)
	dates <- curr.burrow.dat$date 
	
	#split into morning and evening data
	gps.dat.morn <- gps[which(hour(gps$ts) <= morning.evening.thresh),]
	gps.dat.eve <- gps[which(hour(gps$ts) > morning.evening.thresh),]

	if(!is.null(group)){
		curr.gps.dat.morn <- gps.dat.morn[which((as.Date(gps.dat.morn$VisitDate) %in% (dates+1)) & (gps.dat.morn$group==group)),]
		curr.gps.dat.eve <- gps.dat.eve[which((as.Date(gps.dat.eve$VisitDate) %in% (dates)) & (gps.dat.eve$group==group)),]
	} else{
		curr.gps.dat.morn <- gps.dat.morn[which((as.Date(gps.dat.morn$VisitDate) %in% (dates+1))),]
		curr.gps.dat.eve <- gps.dat.eve[which((as.Date(gps.dat.eve$VisitDate) %in% (dates))),]
	}

	#Get plot area
	minlat <- min(min(curr.gps.dat.morn$Latitude,na.rm=T),min(curr.gps.dat.eve$Latitude,na.rm=T))
	minlon <- min(min(curr.gps.dat.morn$Longitude,na.rm=T),min(curr.gps.dat.eve$Longitude,na.rm=T))
	maxlat <- max(max(curr.gps.dat.morn$Latitude,na.rm=T),max(curr.gps.dat.eve$Latitude,na.rm=T))
	maxlon <- max(max(curr.gps.dat.morn$Longitude,na.rm=T),max(curr.gps.dat.eve$Longitude,na.rm=T))
	e <- extent(c(minlon,maxlon,minlat,maxlat))
	g <- gmap(e, type='satellite',lonlat=T)

	#plot the data from each day as lines
	plot(g)
	if(morn){
		
		for(i in 1:length(dates)){
			data <- curr.gps.dat.morn[which(as.Date(curr.gps.dat.morn$ts)==(dates[i]+1)),]
			if(nrow(data)>0){
				times <- hour(data$ts)*60 + minute(data$ts)
				ranks <- rank(times)
				data <- data[ranks,]
				lines(data$Longitude,data$Latitude,col=morning.cols[times],lwd=1)
				points(data$Longitude,data$Latitude,col=morning.cols[times],pch=19,cex=0.5)
			}
		}
		#plot the burrow location
		points(curr.burrow.dat$lon[1],curr.burrow.dat$lat[1],pch=19,col=burrow.col,cex=1)
	}
	
	if(eve){
		for(i in 1:length(dates)){
			data <- curr.gps.dat.eve[which(as.Date(curr.gps.dat.eve$ts)==dates[i]),]
			if(nrow(data)>0){
				times <- hour(data$ts)*60 + minute(data$ts)
				ranks <- rank(times)
				data <- data[ranks,]
				lines(data$Longitude,data$Latitude,col=evening.cols[times],lwd=1)
				points(data$Longitude,data$Latitude,col=evening.cols[times],pch=19,cex=0.5)
			}
		}
		#plot the burrow location
		points(curr.burrow.dat$lon[1],curr.burrow.dat$lat[1],pch=19,col=burrow.col,cex=1)
	}	
}


#Function to make a heatmap of data around a burrow
#Inputs:
#	burrow: the burrow ID number
#	group: the group ID number (or NULL, if plotting data from all groups)
#	time.of.day: string specifying 'morn','eve', or 'both' (which time(s) of day to output - defaults to both)
#	burrow.use.over.time: data frame of burrow data over time
#	gps: data frame of GPS data
#	grid.size: size of the raster grid to use
#	morning.evening.thresh: dividing hour between morning and evening (defaults to 13, meaning 13:00)
#	verbose: whether to print function progress
#Outpus:
#	out: a list containing the following elements
#		morn.rast, eve.rast: rasters for morning and evening data (when group was at the specified burrow) with number representing frequency of occurrence in that spot
#		group, burrow: as in inputs
#		dates: list of dates that group used that burrow (corresponding to evenings) 
#		morn.dates, eve.dates: dates from which trajectories were used for both morning and evening data
burrow.heatmap <- function(burrow,group,group.diam=20,include.lines=F,line.width=10,time.of.day='both',burrow.use.over.time=burrow.use.over.time,gps=gps,grid.size=5,morning.evening.thresh=13,verbose=F){
	
	if(verbose){
		print(paste('creating burrow heatmap for time.of.day =',time.of.day,', burrow =',burrow,', group=',group))
	}
	
	#whether to plot morning or evening data or both (default)
	morn <- eve <- T
	if(time.of.day %in% c('morning','morn','Morn','Morning')){
		eve <- F
	}
	if(time.of.day %in% c('evening','eve','Eve','Evening')){
		morn <- F
	}
	
	#get data for the specified burrow (and group, if needed)
	if(!is.null(group)){
		curr.burrow.dat <- burrow.use.over.time[which(burrow.use.over.time$burrow==burrow & burrow.use.over.time$group==group),] 
	} else{
		curr.burrow.dat <- burrow.use.over.time[which(burrow.use.over.time$burrow==burrow),]
	}

	#get dates that burrow was used (evenings)
	dates <- curr.burrow.dat$date 
	
	#split into morning and evening data
	gps.dat.morn <- gps[which(hour(gps$ts) <= morning.evening.thresh),]
	gps.dat.eve <- gps[which(hour(gps$ts) > morning.evening.thresh),]

	if(!is.null(group)){
		curr.gps.dat.morn <- gps.dat.morn[which((as.Date(gps.dat.morn$VisitDate) %in% (dates+1)) & (gps.dat.morn$group==group)),]
		curr.gps.dat.eve <- gps.dat.eve[which((as.Date(gps.dat.eve$VisitDate) %in% (dates)) & (gps.dat.eve$group==group)),]
	} else{
		curr.gps.dat.morn <- gps.dat.morn[which((as.Date(gps.dat.morn$VisitDate) %in% (dates+1))),]
		curr.gps.dat.eve <- gps.dat.eve[which((as.Date(gps.dat.eve$VisitDate) %in% (dates))),]
	}
	
	#get appropriate data
	if(eve & morn){
		curr.gps <- rbind(curr.gps.dat.morn,curr.gps.dat.eve)
	} else{
		if(eve){
			curr.gps <- curr.gps.dat.eve
		}
		if(morn){
			curr.gps <- curr.gps.dat.morn
		}
	}
	
	#get minimum and maximum extent of data
	minx <- floor(min(curr.gps$Easting,na.rm=T))
	miny <- floor(min(curr.gps$Northing,na.rm=T))
	maxx <- ceiling(max(curr.gps$Easting,na.rm=T))
	maxy <- ceiling(max(curr.gps$Northing,na.rm=T))
	
	#get an image
	morn.rast <- eve.rast <- raster(ext=extent(minx,maxx,miny,maxy),vals=0,res=grid.size)
	eve.dates <- morn.dates <- c()
	for(i in 1:length(dates)){
		
		if(verbose && (i %% 100 == 0)){
			print(paste(i,'/',length(dates)))
		}
		if(morn){
			data <- curr.gps.dat.morn[which(as.Date(curr.gps.dat.morn$ts)==(dates[i])+1),]
			if(nrow(data)>0){
				morn.dates <- c(morn.dates,dates[i]+1)
				morn.rast <- morn.rast + trajectory.heatmap(data$Easting,data$Northing,minx,maxx,miny,maxy,grid.size=grid.size,diam=group.diam,include.lines=include.lines,line.width=line.width)
			}
		}
		if(eve){
			data <- curr.gps.dat.eve[which(as.Date(curr.gps.dat.eve$ts)==dates[i]),]
			if(nrow(data)>0){
				eve.dates <- c(eve.dates,dates[i])
				eve.rast <- eve.rast + trajectory.heatmap(data$Easting,data$Northing,minx,maxx,miny,maxy,grid.size=grid.size,diam=group.diam,include.lines=include.lines,line.width=line.width)
			}
		}

	}
	out <- list(morn.rast=morn.rast,eve.rast=eve.rast,group=group,burrow=burrow,dates=dates,morn.dates=morn.dates,eve.dates=eve.dates)
	return(out)
}

#Make a heatmap (raster) of the points in a trajectory (and optionally the lines connecting them)
#Function is called by burrow.heatmap
#Inputs:
#	xs, ys: x and y coordinates of trajectory (vectors)
#	minx, maxx, miny, maxy: extent of the raster
#	grid.size: pixel size for the raster
#	diam: diamaeter of the group to use (how big of a 'blob' to use in the heatmap)
#	include.lines: whether to include lines between points in the trajectory
#	line.width: width of lines to use, if included
#Outputs:
#	rast: a raster containing the heatmap (1's if the trajectory crosses that space, 0 if not)
trajectory.heatmap <- function(xs,ys,minx,maxx,miny,maxy,grid.size=5,diam=50,include.lines=F,line.width=10){
	img <- raster(extent(minx,maxx,miny,maxy),res=grid.size)
	rast <-distanceFromPoints(img,cbind(xs,ys)) < (diam/2)
	if(include.lines){
		line.rast <- rasterize(spLines(cbind(xs,ys)),img)
		dist.lines <- distance(line.rast) < line.width
		rast <- (rast + dist.lines) > 0
	}
	rast[rast==FALSE] <- 0
	rast[rast==TRUE] <- 1
	return(rast)	
}

#Function to convert data in the form xx:xx:xx.xxx OR xx:xx.xxx to seconds (uses lubridate)
#Inputs:
#	instr: input string
#Outputs:
#	secs: number of seconds (numeric)
to.secs <- function(instr){
	if(length(grep(':..:',instr))>0){
		secs <- hms(instr)
	} else{
		secs <- ms(instr)
	}
	secs <- seconds(secs)
	return(secs)
}


library(dismo)

#Generate a bunch of images (PNGs) showing trajectory data over a given time period (optionally on top of a Google Earth map)
#These can then be ffmpeg'ed together into a video.
#Inputs:
#	lats: [N x T matrix] of latitude coordinates (N individuals, T timesteps)
#	lons: [N x T matrix] of longitude coordinates (N individuals, T timesteps)
#	start.time: [numeric] time index at which to start the video
#	end.time: [numeric] time index at which to end the video
#	tail.time: [numeric] number of previous time steps to plot as a "tail" which trails the point showing the current location
#	base.dir: [string] directory in which to store the folder of outputted images
#	colors: OPTIONAL [vector of strings of length N] to specify different colors for different individuals (if NULL, default is to color individuals randomly using a rainbow palette)
#	on.map: OPTIONAL [boolean]. If TRUE, trajectories will be plotted on a Google Earth map. If false, they will be plotted on a white background.
#	ind.names: [vector of strings] giving names of individuals
#	plot.legend: [boolean] whether to show a legend
#	times: [vector of times] to show on map
#	show.all.past: [boolean] whether to show all the past data in a lighter color 
#	calls: OPTIONAL [data frame] of calls to plot on top of trajectories (defaults to NULL)
#	call.persist.time: OPTIONAL [numeric] how long the calls are shown for (by default, this matches tail time)
#	call.types: OPTIONAL [data frame] of call types and colors / symbols to plot them in (defaults to black circles for all)
#	zoom: OPTIONAL [numeric] from 8 to 21 (21 is fully zoomed in)
# places: either NULL (if no landmark places) or a data frame containing:
#   $name: name of each place
#   $lon, lat: lon and lat coordinates of each place
#   $color: color to plot in
#   $cex: size of the point
#   $pch: pch of the point
# traces: either NULL (if no 'traces' i.e. paths) or a list of n x 2 matrices containing lon/lat coordinates of lines  
#Outputs:
#	a folder of PNG images inside base.dir
trajectories.movie <-function(lats,lons,start.time,end.time,step=1,tail.time=9,base.dir='/Users/astrandb/Desktop',colors=NULL,on.map=TRUE,ind.names=NULL,plot.legend=F,times=NULL,show.all.past=FALSE,show.scale.bar=TRUE,scale.bar.len=1000,utm.zone=36,southern_hemisphere=T,scale.bar.text=NULL,scale.bar.text.offset=NULL,playback.time=NULL,playback.lat=NULL,playback.lon=NULL,inds=NULL,calls=NULL,call.persist.time=NULL,call.types=NULL,zoom=20,past.step=1,places=NULL, traces = NULL){
	
	#create directory in which to store images
	dir.create(paste(base.dir,'/seq',start.time,'-',end.time,sep=''))
	
	#number of individuals
	N <- dim(lats)[1]
	
	#get colors if not specified
	if(is.null(colors)){
		colors <- rainbow(N)
	}
	
	#get lats and lons to use
	if((start.time - tail.time) > 1){
		currlats = lats[,(start.time-tail.time):end.time]
		currlons = lons[,(start.time-tail.time):end.time]
	} else{
		currlats = lats[,1:end.time]
		currlons = lons[,1:end.time]
	}
	
	
	#get times vector
	if(!is.null(times)){
		if((start.time - tail.time) > 1){
			currtimes <- times[(start.time-tail.time):end.time]
		} else{
			currtimes <- times[1:end.time]
		}
	}
	
	if((start.time - tail.time) < 1){
		min.time <- 1
	} else{
		min.time <- start.time - tail.time
	
	}
	#get map boundaries - depending on your scale you may need to adjust the .0005's
	
	if(is.null(inds)){
		ind.idxs <- seq(1,N)
	} else{
		ind.idxs <- inds
	}
	
	xmin = quantile(lons[ind.idxs,min.time:end.time],0.0001,na.rm=T)
	xmax = quantile(lons[ind.idxs,min.time:end.time],0.9999,na.rm=T)
	ymin = quantile(lats[ind.idxs,min.time:end.time],0.0001,na.rm=T)
	ymax = quantile(lats[ind.idxs,min.time:end.time],0.9999,na.rm=T)
	
	#set map extent
	e<-extent(xmin,xmax,ymin,ymax)
	
	#get map if needed
	if(on.map){
		g<-gmap(e,type='satellite',lonlat=TRUE)
	}
	
	#get legend info if needed
	if(plot.legend){
		legend.x <- e@xmin + (e@xmax-e@xmin)/100
		legend.y <- e@ymax - (e@ymax-e@ymin)/100
	}
	
	#get time colors
	time.cols.palette <- colorRampPalette(c('black','black','orange','yellow','yellow','yellow','orange','black','black'))
	time.cols <- time.cols.palette(24)
	
	#if scale bar needed, source utm to lat lon, and get easts and norths
	if(show.scale.bar){
		eastsNorths <- latlon.to.utm(cbind(xmax,ymin),southern_hemisphere=southern_hemisphere,utm.zone=utm.zone)
		lonsLats <- utm.to.latlon(cbind(eastsNorths[,1]-scale.bar.len,eastsNorths[,2]),southern_hemisphere=southern_hemisphere,utm.zone=utm.zone)
		scale.minlon <- lonsLats[,1]
		scale.minlat <- ymin
		scale.maxlon <- xmax
		scale.maxlat <- ymin
		text.lon <- (scale.minlon + scale.maxlon)/2
		text.lat <- scale.minlat + scale.bar.text.offset
		lonsLats2 <- utm.to.latlon(cbind(eastsNorths[,1]-scale.bar.len/2,eastsNorths[,2]-scale.bar.text.offset),southern_hemisphere=southern_hemisphere,utm.zone=utm.zone)
	}
	
	#idx<-1
	idx <- 1
	prev.lons <- lons[,start.time]
	prev.lats <- lats[,start.time]
	for(i in seq(start.time,end.time,step)){

		#get lats and lons for past and present data
		x<-lons[,i]
		y<-lats[,i]
		
		for(j in 1:N){
			if(!is.na(x[j])){
				prev.lons[j] <- x[j]
				prev.lats[j] <- y[j]
			}
		}
		
		#get data for tail
		if(tail.time > 0){
			if((i - tail.time) < 1){
				past.idxs <- 1:i
			} else{
				past.idxs <- (i-tail.time):i
			}
			x.past<-as.matrix(lons[,past.idxs])
			y.past<-as.matrix(lats[,past.idxs])
		}
		
		#get data for translucent tail (all past data)
		if(show.all.past){
			x.past.all <- as.matrix(lons[,seq(1,i,past.step)])
			y.past.all <- as.matrix(lats[,seq(1,i,past.step)])
		}
		
		#make figure
		filename = paste(base.dir,'/seq',start.time,'-',end.time,'/',idx,'.png',sep='')
		png(file=filename,width=10,height=4,units='in',res=300)
		par(mar=c(0,0,0,0))
		
		#initialized plot or map
		if(on.map){
			plot(g,interpolate=TRUE)
		} else{
			par(bg='black')
			plot(NULL,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxt='n',yaxt='n',xlab='',ylab='',bg='black',asp=1)
		}
		
		#plot all past data (if required)
		if(show.all.past){
			for(j in 1:N){
				points(x.past.all[j,],y.past.all[j,],col=adjustcolor(colors[j],alpha.f=0.1),bg=adjustcolor(colors[j],alpha.f=0.1),pch=19,cex=0.05)
			}
			
		}
		
		#plot "tails" (past locations)
		for(j in 1:N){
			lines(x.past[j,],y.past[j,],col=colors[j],lwd=2)
			#points(x.past[j,],y.past[j,],col=colors[j],cex=0.7,pch=19,bg=colors[j])
		}
		
		#add a playback (if needed)
		if(!is.null(playback.time)){
			if(abs(i-playback.time) < 1200){
				if(i < playback.time){
					points(playback.lon,playback.lat,pch=15,cex=2,col='black')
				}
				if(i >= playback.time){
					if(abs(i-playback.time)<10){
						points(playback.lon,playback.lat,pch=15,cex=2,col='yellow')
						points(playback.lon,playback.lat,pch=8,cex=2,col='black')
					} else{
						points(playback.lon,playback.lat,pch=15,cex=2,col='yellow')
					}
				}
			}		
		}
		
		#plot a legend
		if(plot.legend){
			if(!is.null(ind.names) & is.null(call.types)){
				legend(legend.x,legend.y,as.character(ind.names),col=colors,fill=colors,border='white',text.col='white',box.col='white',bty='n')
			}
			if(!is.null(ind.names) & !is.null(call.types)){
				legend('bottomleft',c(as.character(ind.names),'',as.character(call.types$type)),col=c(colors,'black',as.character(call.types$col)),pch=c(rep(19,length(ind.names)+1),call.types$sym),border='white',text.col='white',box.col='white',bty='n')
			}
		}
		
		#plot the time
		if(!is.null(times)){
			#text(x=xmin,y=ymin,labels=times[i],col=time.cols[hour(times[i])+1],cex=1)
			text(x=xmin+.004,y=ymax-.0001,labels=times[i],col='white',cex=1)
		}
		
		#make a scale bar
		if(show.scale.bar){
			#lines(c(scale.minlon,scale.maxlon),c(scale.minlat,scale.maxlat),lwd=2,col='white')
			lines(c(scale.minlon,scale.maxlon),c(scale.minlat+.00002,scale.maxlat+.00002),lwd=2,col='white')
			#text(x=text.lon,y=text.lat,labels=scale.bar.text,col='white')
			text(x=text.lon,y=scale.minlat+.00001,labels=scale.bar.text,col='white')
		}
		
		#plot current locations
		#for(j in 1:N){
		#	points(prev.lons[j],prev.lats[j],col=adjustcolor(colors[j],alpha.f=0.2),bg=adjustcolor(colors[j],alpha.f=0.2),pch=19,cex=1)
		#}
		if(on.map){
			points(x,y,pch=21,cex=1.5,col='white',bg=colors)
		} else{
			points(x,y,pch=21,cex=.5,col=colors,bg=colors)
		}
		
		#plot calls in the past
		if(!is.null(calls)){
  		calls.curr <- calls[which(is.na(calls$nonfoc) & calls$t.idx < i & calls$t.idx >= (i-call.persist.time)),]
  		if(nrow(calls.curr)>0){
  			calls.lon <- lons[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
  			calls.lat <- lats[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
  			calls.type.idxs <- match(calls.curr$call.type,call.types$type)
  			calls.col = as.character(call.types$col[calls.type.idxs])
  			calls.sym = call.types$sym[calls.type.idxs]
  			points(calls.lon,calls.lat,col=calls.col,pch=calls.sym,cex=1)
  		}
		}
		
		#plot current calls
		if(!is.null(calls)){
			calls.curr <- calls[which(calls$t.idx == i & is.na(calls$nonfoc)),]
			if(nrow(calls.curr)>0){
				calls.lon <- lons[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
				calls.lat <- lats[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
				calls.type.idxs <- match(calls.curr$call.type,call.types$type)
				calls.col = as.character(call.types$col[calls.type.idxs])
				calls.sym = call.types$sym[calls.type.idxs]
				points(calls.lon,calls.lat,col=calls.col,pch=calls.sym,cex=2)
			}
		}
		
		if(!is.null(traces)){
		  for(tr in 1:length(traces)){
		    xtr <- traces[[tr]][,1]
		    ytr <- traces[[tr]][,2]
		    lines(xtr,ytr,lwd=2,col='white')
		  }
		}
		
		if(!is.null(places)){
		  points(places$lon,places$lat,col=places$col,pch=places$pch,cex=places$cex,lwd=places$lwd)
		}
		
		dev.off()
		
		idx<-idx+1

	}
	
}

#Generate a circular histogram (or a heat map)
#INPUTS:
# rs: [numeric] vector of radii to use 
# thetas: [numeric] vector of angles to use
# vals: [numeric] vector of values to go along with rs and thetas (or xs and ys) that will be summarized over the heat map using summary.stat (default NULL gives a regular histogram)
# xs: [numeric] vector of x coordinates to use (optional, will only be used if rs and thetas are NULL)
# ys: [numeric] vector of y coordinates to use (optional, will only be used if rs and thetas are NULL)
# r.bins: [numeric] vector of radius bins to use
# theta.bins: [numeric] vector of radius bins to use
# summary.stat: [chr] specification of either 'mean' or 'median' for how to summarize the values
# cols: vector of hex color values
# zlim: [numeric] vector of length 2 giving minimum and maximum values to scale colors by 
heatmap.circ <- function(rs = NULL, thetas = NULL, vals = NULL, xs = NULL, ys = NULL, r.bins = NULL, theta.bins = NULL, summary.stat = 'mean',cols=NULL,zlim=NULL,xlab='',ylab=''){
  
  #convert to polar coordinates if needed
  if(is.null(rs) | is.null(thetas)){
    
    #check to make sure xs and ys have been input, otherwise throw an error
    if(is.null(ys) | is.null(xs)){
      error('please specify either rs and thetas or xs and ys')
    }
    
    rs <- sqrt(xs^2 + ys^2)
    thetas <- atan2(ys,xs)

  }
  
  #default values of bins
  if(is.null(theta.bins)){
    theta.bins <- seq(-pi,pi,pi/8)
  }
  if(is.null(r.bins)){
    r.bins <- quantile(rs,seq(0,.95,.05),na.rm=T)
  }
  
  #get colors (default to viridis or red/blue depending on whether vals are included or not)
  if(is.null(cols)){
    if(is.null(vals)){
      cols <- viridis(1024)
    } else{
      col.palette <- colorRampPalette(c('red','white','blue'))
      cols <- col.palette(1024)
    }
  }
  
  #compute z based on probabilities of being in each bin or a summary statistic of 'vals' within each bin
  out.mat <- ns <- matrix(NA,nrow=length(r.bins)-1,ncol=length(theta.bins)-1)
  for(i in 1:(length(r.bins)-1)){
    for(j in 1:(length(theta.bins)-1)){
      idxs <- which((rs >= r.bins[i]) & (rs < r.bins[i+1]) & (thetas >= theta.bins[j]) & (thetas < theta.bins[j+1]))
      if(is.null(vals)){
        ns[i,j] <- length(idxs)
        out.mat[i,j] <- length(idxs)
      } else{
        if(summary.stat=='mean'){
          ns[i,j] <- length(idxs)
          out.mat[i,j] <- mean(vals[idxs],na.rm=T)
        } else{
          if(summary.stat =='median'){
            ns[i,j] <- length(idxs)
            out.mat[i,j] <- median(vals[idxs],na.rm=T)
          } else{
            error('please enter a valid summary stat (mean or median currently supported)')
          }
        }
      }
    }
  }
  
  print(out.mat)
  
  #normalize histogram (if it's a regular histogram)
  if(is.null(vals)){
    out.mat <- out.mat / sum(out.mat,na.rm=T)
  }
  
  if(is.null(zlim)){
    if(is.null(vals)){
      zlim = c(0,max(out.mat,na.rm=T))
    } else{
      zlim = c(min(out.mat,na.rm=T),max(out.mat,na.rm=T))
    }
    
  }
  
  #make the plot
  #quartz()
  xlim <- c(-max(r.bins),max(r.bins))
  ylim <- c(-max(r.bins),max(r.bins))
  plot(NULL,xlim=xlim,ylim=ylim,asp=1,xlab=xlab,ylab=ylab,bty='n',xaxt='n',yaxt='n')
  for(i in 1:(length(r.bins)-1)){
    for(j in 1:(length(theta.bins)-1)){
      dat <- out.mat[i,j]
      col <- cols[(dat - zlim[1])/(zlim[2] - zlim[1])*(length(cols)-1)+1]
      r0 <- r.bins[i]
      rf <- r.bins[i+1]
      ang0 <- theta.bins[j]
      angf <- theta.bins[j+1]
      d.angs <- seq(ang0,angf,(angf - ang0)/100)
      x1 <- r0*cos(ang0)
      x2 <- rf*cos(d.angs)
      x3 <- r0*cos(angf)
      x4 <- r0*cos(rev(d.angs))
      y1 <- r0*sin(ang0)
      y2 <- rf*sin(d.angs)
      y3 <- r0*sin(angf)
      y4 <- r0*sin(rev(d.angs))
      polygon(c(x1,x2,x3,x4),c(y1,y2,y3,y4),col=col,border=NA)
    }
  }
  
  #colorbar.plot(x=0,y=-max(r.bins)*1.2,strip=seq(zlim[1],zlim[2],length.out = length(cols)),col=cols,strip.width = 0.05,strip.length = strip.len,horizontal = T,)
  image.plot(legend.only=T,zlim=zlim,col=cols,horizontal=T)
  
  lines(c(0,max(r.bins,na.rm=T)),c(0,0),col='black',lwd=2)
  text(as.character(round(max(r.bins,na.rm=T),digits = 1)),x=max(r.bins,na.rm=T),y=0,pos=4)
  return(list(out.mat,ns))
}


#FUNCTION: Gets the call rate for a focal meerkat (caller) at time t over the past time.win seconds
#INPUTS:
# t: time of computation (index)
# time.win : time window into the past 
# caller: index of focal meerkat
# calls: data frame of call data
# times: time strings 
#OUTPUTS:
# n.calls.per.minte: the number of calls per minute
get.call.rate <- function(t,time.win,caller,calls,times,call.type='CC'){
  #get start and end time of window over which calls will be found
  start.time <- t-time.win
  end.time <- t
  
  #check whether we are in a range where there is data, otherwise return NA
  
  #get the data of the current time point
  curr.date <- date(t)
  start.ind.date.row <- which(calls$date==curr.date & calls$ind.idx==caller & calls$call.name=='START')
  end.ind.date.row <- which(calls$date==curr.date & calls$ind.idx==caller & calls$call.name=='END')
  
  #if there is no data for that individual on that date, return NA
  if(length(start.ind.date.row)==0){
    return(NA)
  }
  
  #otherwise, we need to check whether our time is in the right range
  t.label.begin <- calls$t0.gps[start.ind.date.row]
  t.label.end <- calls$t0.gps[end.ind.date.row]
  
  if(length(t.label.begin)>0){
    current.times <- c()
    for(i in 1:length(t.label.begin)){
      current.times <- round(as.numeric(c(current.times,t.label.begin[i]:t.label.end[i])))
    }
  }
  
  if((round(as.numeric(t)) %in% current.times & round(as.numeric(t-time.win)) %in% current.times)){
    
    #find rows within time window associated with the caller
    rows.in.window <- which(calls$t0.gps <= end.time & calls$t0.gps > start.time & calls$ind.idx==caller)
    
    #compute number of calls in those rows then convert to number per minute
    n.calls.per.minute <- sum(calls$call.type[rows.in.window]==call.type)/time.win*60
    
    #return the output
    return(n.calls.per.minute)
  } else{
    return(NA)
  }
  
}

#Function to get call rate over a time window t.win ending at t.idx
#This function assumes that user has already checked that the region t.idx-t.win to t.idx
#is in the labeled period, and that call.time.idxs correspond to the call type you want
get.call.rate.fast <- function(call.time.idxs,t.idx,t.win){
  tot.calls <- sum((call.time.idxs >= (t.idx-t.win)) & (call.time.idxs < t.idx),na.rm=T)
  return(tot.calls / t.win * 60)
}


#Function to determine angle of toward/away based on the time when focal leaves a circle of radius R
#INPUTS:
# xs, ys: [N x T matrix] of x and y coords of all individuals
# i: focal individual (the "mover")
# j: nonfocal individual (the "caller")
# t0: initial time index
# R: distance they have to move before getting the angle (in meters)
#OUTPUTS:
# ang: angle towards or away (in degrees)
toward.away.angle.spatial <- function(xs,ys,i,j,t0,R){
  
  #get dimensions of matrices
  n.inds <- dim(xs)[1]
  n.times <- dim(xs)[2]
  
  #get x and y coords for each individual at all times after t0
  xi <- xs[i,t0:n.times]
  yi <- ys[i,t0:n.times]
  xj <- xs[j,t0:n.times]
  yj <- ys[j,t0:n.times]
  
  #get initial vector pointing from i to j
  dx0 <- xi[1] - xj[1]
  dy0 <- yi[1] - yj[1]
  
  #get initial position of i
  xi0 <- xi[1]
  yi0 <- yi[1]
  
  #for each time in the future, get distance from original point, stop when exceeds R
  curr.dist <- 0
  t.idx <- 1
  while(curr.dist < R | (t.idx == length(xi))){
    curr.dist <- sqrt((xi0 - xi[t.idx])^2 + (yi0 - yi[t.idx])^2)
    t.idx <- t.idx + 1
  }
  
  #position at final time
  xif <- xi[t.idx]
  yif <- yi[t.idx]
  
  #vector pointing to final position from initial position
  dxi <- xif - xi0
  dyi <- yif - yi0
  
  #get angle between vectors via law of cosines
  num <- dx0 * dxi + dy0 * dyi #numerator for law of cosines dot product
  len.d0 <- sqrt( dx0^2 + dy0^2 ) #length of initial vector between meerkats
  len.di <- sqrt( dxi^2 + dyi^2 ) #length of the vector of movement of the focal i
  cosang <- num / (len.d0 * len.di) #cosine of angle
  ang <- acos(cosang)*180/pi
  
  return(ang)
  
}

#Get the index to the minimum in a vector with ties broken at random, ignoring NAs
which.min.random <- function(x){
  
  #find minimum value
  min.val <- min(x,na.rm=T)
  
  if(is.infinite(min.val)){
    return(NA)
  }
  
  whiches <- which(x==min.val)
  
  if(length(whiches)>1){
    return(sample(whiches,1))
  } else{
    return(whiches[1])
  }
  
}


#UPDATED VERSION OF GMAP
gmap <- function (x, exp = 1, type = "terrain", filename = "", style = NULL, 
          scale = 1, zoom = NULL, size = c(640, 640), rgb = FALSE, 
          lonlat = FALSE, ...) 
{
  if (!requireNamespace("rgdal")) {
    stop("rgdal not available")
  }
  if (!type %in% c("roadmap", "satellite", "hybrid", "terrain")) {
    warning("type should be: roadmap, satellite, hybrid, or terrain: Terrain chosen by default")
    type <- "terrain"
  }
  mxzoom <- function(latrange, lonrange, size = size) {
    SinPhi = sin(latrange * pi/180)
    normX = lonrange/180
    normY = (0.5 * log(abs((1 + SinPhi)/(1 - SinPhi))))/pi
    MaxZoom.lon <- floor(1 + log2(abs(size[1]/256/diff(normX))))
    MaxZoom.lat <- floor(1 + log2(abs(size[2]/256/diff(normY))))
    return(c(MaxZoom.lat = MaxZoom.lat, MaxZoom.lon = MaxZoom.lon))
  }
  ll2XY <- function(lat, lon, zoom) {
    SinPhi = sin(lat * pi/180)
    normX = lon/180
    normY = (0.5 * log((1 + SinPhi)/(1 - SinPhi)))/pi
    Y = (2^zoom) * ((1 - normY)/2)
    X = (2^zoom) * ((normX + 1)/2)
    x = 256 * (X - floor(X))
    y = 256 * (Y - floor(Y))
    return(list(Tile = cbind(X = floor(X), Y = floor(Y)), 
                Coords = cbind(x = x, y = y)))
  }
  xy2ll <- function(MyMap, X, Y) {
    lat.center <- MyMap[[1]]
    lon.center <- MyMap[[2]]
    zoom <- MyMap[[3]]
    mycenter <- ll2XY(lat.center, lon.center, zoom)
    x <- mycenter$Tile[, "X"] + (X + mycenter$Coords[, "x"])/256
    y <- mycenter$Tile[, "Y"] - (Y - mycenter$Coords[, "y"])/256
    ytilde <- 1 - y/2^(zoom - 1)
    yy = (exp(2 * pi * ytilde) - 1)/(exp(2 * pi * ytilde) + 
                                       1)
    ShiftLat <- function(yy) {
      n = c(-1, 0, 1)
      lat = 2 * pi * (n) + asin(yy)
      lat <- lat[which(lat <= pi/2 & lat > -pi/2)]
      lat <- 180 * lat/pi
      return(lat)
    }
    lat <- sapply(yy, ShiftLat)
    lon = 180 * (x/2^(zoom - 1) - 1)
    return(cbind(lat = lat, lon = lon))
  }
  tile2r <- function(points, center) {
    X <- 256 * (points$Tile[, "X"] - center$Tile[, "X"]) + 
      (points$Coords[, "x"] - center$Coords[, "x"])
    Y <- -256 * (points$Tile[, "Y"] - center$Tile[, "Y"]) - 
      (points$Coords[, "y"] - center$Coords[, "y"])
    return(list(X = X, Y = Y))
  }
  #gurl <- "http://maps.googleapis.com/maps/api/staticmap?"
  gurl <- 'https://maps.googleapis.com/maps/api/staticmap?'
  if (is.character(x)) {
    x <- geocode(x, oneRecord = TRUE)
    if (is.na(x$latitude)) {
      stop("location not found")
    }
    x <- extent(as.vector(as.matrix(x[5:8])))
  }
  else {
    prj <- projection(x, asText = TRUE)
    if (isLonLat(prj)) {
      x <- extent(x)
    }
    else {
      if (is.na(prj)) {
        bb <- extent(x)
        extLL <- (bb@xmin > -366 & bb@xmax < 366 & bb@ymin > 
                    -90.1 & bb@ymax < 90.1)
        if (extLL) {
          x <- bb
        }
        else {
          rad <- 6378137
          p <- t(bbox(x))
          p[, 2] <- pi/2 - 2 * atan(exp(-p[, 2]/rad))
          p[, 1] <- p[, 1]/rad
          p <- p/(pi/180)
          x <- extent(p[1, 1], p[2, 1], p[1, 2], p[2, 
                                                   2])
        }
      }
      else {
        x <- extent(projectExtent(x, "+proj=longlat +datum=WGS84"))
      }
    }
  }
  e <- x * exp
  e@xmin <- max(-180, e@xmin)
  e@xmax <- min(180, e@xmax)
  e@ymax <- min(89, e@ymax)
  e@ymin <- max(-89, e@ymin)
  lonR <- c(e@xmin, e@xmax)
  latR <- c(e@ymin, e@ymax)
  if (is.null(zoom)) {
    zoom <- min(mxzoom(latR, lonR, size))
  }
  center <- c(mean(latR), mean(lonR))
  ll <- ll2XY(latR[1], lonR[1], zoom)
  ur <- ll2XY(latR[2], lonR[2], zoom)
  cr <- ll2XY(center[1], center[2], zoom)
  ll.Rcoords <- tile2r(ll, cr)
  ur.Rcoords <- tile2r(ur, cr)
  s1 <- 2 * max(c(ceiling(abs(ll.Rcoords$X)), ceiling(abs(ur.Rcoords$X)))) + 
    1
  s2 <- 2 * max(c(ceiling(abs(ll.Rcoords$Y)), ceiling(abs(ur.Rcoords$Y)))) + 
    1
  if (s1 > s2) {
    size[2] <- as.integer(size[2] * s2/s1)
  }
  else {
    size[1] <- as.integer(size[1] * s1/s2)
  }
  s <- paste(size, collapse = "x")
  ctr <- paste(center, collapse = ",")
  gurl <- paste(gurl, "center=", ctr, "&zoom=", zoom, "&size=", s, "&maptype=", type,"&key=", key,"&format=gif&sensor=false&scale=", scale, sep = "")
  if (!is.null(style)) {
    style <- gsub("\\|", "%7C", style)
    style <- gsub(" ", "", style)
    gurl <- paste(gurl, "&style=", style, sep = "")
  }
  filename <- trim(filename)
  if (filename == "") {
    filename <- rasterTmpFile()
  }
  extension(filename) <- "gif"
  download.file(gurl, filename, mode = "wb", quiet = TRUE)
  MyMap <- list(lat.center = center[1], lon.center = center[2], 
                zoom = zoom)
  bb <- list(ll = xy2ll(MyMap, X = -size[1]/2 + 0.5, Y = -size[2]/2 - 
                          0.5), ur = xy2ll(MyMap, X = size[1]/2 + 0.5, Y = size[2]/2 - 
                                             0.5))
  r <- raster(filename, warn = FALSE)
  ext <- extent(bb$ll[2], bb$ur[2], bb$ll[1], bb$ur[1])
  p <- t(bbox(raster(ext))) * pi/180
  rad <- 6378137
  p[, 2] <- log(tan(p[, 2]) + (1/cos(p[, 2])))
  p <- p * rad
  extent(r) <- extent(as.vector(p))
  projection(r) <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
  try(hdr(r, format = "worldfile", extension = ".gfw"))
  extension(filename) <- "prj"
  rgdal::showWKT(projection(r), file = filename, morphToESRI = TRUE)
  if (lonlat) {
    ct <- r@legend@colortable
    r <- projectRaster(r, crs = "+proj=longlat +datum=WGS84", 
                       method = "ngb")
    r@legend@colortable <- ct
  }
  if (rgb) {
    d <- t(col2rgb(r@legend@colortable))
    d <- data.frame(id = 0:255, d)
    r <- subs(r, d, which = 2:4)
  }
  return(r)
}

#Computes turning leadership score 
#This measure consists of these steps:
#1. Compute the projection of the unit vector pointing from individual f (follower) at time t0 to individual l (leader) at time t0 along the direction perpendicular to the group direction of motion (call this proj_fl)
#2. Compute the projection of the unit vector pointing from individual f (follower) at time t0 to individual f (follower) at time t0 + dt along the direction perpendicular to the group direction of motion (call this proj_ff)
#3. If leader and follower are on the same side relative to the group centroid and heading, and leader is farther out than follower, then define the 'turning influence' as proj_ff * sign(prof_fl)
#   If leader and follower are on different sides relative to the group centroid and heading, and/or if leader is closer in than follower, the measure is undefined (NA)
#Option: Instead of defining a fixed dt, one can instead define a fixed radius R that a follower moves before proj_ff is computed. To do this, set dt to NULL in the function and specify R (in m)
#
#INPUTS: 
# xs: [n.inds x n.times matrix] giving x position (UTM easting) of all individuals at each time
# ys: [n.inds x n.times matrix] giving y position (UTM northing) of all individuals at each time
# spat.heads: [n.times x 1] vector giving heading (in radians) of the group
# leader: [numeric] index of the individual to be considered the 'leader' in the position matrices
# follower: [numeric OR char] index of the individual to be considered the 'follower' in the position matrices OR specify 'centroid' to use the group centroid as the follower
# t0: [numeric] initial time at which to compute the measure
# dt: [numeric] time lag to use for computing follower direction of motion (or specify NULL to use a spatial threshold instead of a temporal one)
# R: [numeric] distance the follower has traveled from starting point when heading is computed (only used if dt = NULL)
#
#OUTPUTS: 
# score: [numeric] the turning leadership score for leader relative to follower at time t0

turning.leadership.score <- function(xs, ys, spat.heads, leader, follower, t0, dt, R){
  
  #if starting at the end of the times, return NA
  if(t0 >= ncol(xs)){
    return(NA)
  }
  
  #get group centroid position at t0
  
  #if using the centroid as the follower, compute it without the leader included
  if(follower == 'centroid'){
    centr.x <- mean(xs[-leader,t0],na.rm=T)
    centr.y <- mean(ys[-leader,t0],na.rm=T)
    
    #otherwise compute the centroid using all individauls including the leader and follower
  } else{
    centr.x <- mean(xs[,t0],na.rm=T)
    centr.y <- mean(ys[,t0],na.rm=T)
  }
  
  #get x and y for the follower over time (or the centroid excluding leader if leader == 'centroid')
  if(follower=='centroid'){
    xfs <- colMeans(xs[-leader,],na.rm=T)
    yfs <- colMeans(ys[-leader,],na.rm=T)
  } else{
    xfs <- xs[follower,]
    yfs <- ys[follower,]
  }
  
  #if dt = NULL is specified, get dt dynamically based on when the follower has moved a distance R
  if(is.null(dt)){
    
    #initial position of the follower
    xf0 <- xfs[t0]
    yf0 <- yfs[t0]
    
    #maximum time (to avoid running beyond the length of the vector)
    maxdt <- ncol(xs) - t0
    
    #loop forward in time until the distance between initial position and current position is > R
    flag <- F #flag to indicate whether the individual has moved a radius R (will then get set to true)
    for(dt in 1:maxdt){
      curr.dist.sq <- (xfs[t0+dt] - xf0)^2 + (yfs[t0+dt] - yf0)^2 #current squared distance of the follower from initial position
      
      #if any missing data is encountered, return NA for the measure (missing data is not allowed)
      if(is.na(curr.dist.sq)){
        return(NA)
      }
      
      #if the individual hasm oved at least a distance R, break out of the loop and set flag = T
      if(curr.dist.sq >= R^2){
        flag <- T
        break
      }
    }
    
    #if the radius is not reached (loop exits and flag == F), return NA for the score. Otherwise continue
    if(flag == F){
      return(NA)
    }
  }
  
  #get group heading
  group.heading <- spat.heads[t0]
  head.x <- cos(group.heading)
  head.y <- sin(group.heading)
  
  #get vector 90 degrees rotated from group heading (this performs a 90 degree rotation clockwise)
  head.perp.x <- -head.y
  head.perp.y <- head.x
  
  #get vector from group centroid to leader initial position
  v_centr_l.x <- xs[leader,t0] - centr.x
  v_centr_l.y <- ys[leader,t0] - centr.y
  
  #get vector from centroid to follower initial position
  v_centr_f.x <- xfs[t0] - centr.x
  v_centr_f.y <- yfs[t0] - centr.y
  
  #project both centroid to individual vectors onto direction perpendicular to group heading
  leader_centr_proj <- head.perp.x*v_centr_l.x + head.perp.y*v_centr_l.y
  follower_centr_proj <- head.perp.x*v_centr_f.x + head.perp.y*v_centr_f.y
  
  #if encounter any NAs, return NA for the score
  if(any(is.na(c(leader_centr_proj, follower_centr_proj)))){
    return(NA)
  }
  
  #check whether those two proejctions are the same sign and that leader is farther away from centroid along perpendicular than follower (otherwise return NA)
  if(((leader_centr_proj * follower_centr_proj) < 0) | (abs(leader_centr_proj) <= abs(follower_centr_proj))){
    return(NA)
  }
  
  #get unit vector u_ff from follower current location to its future location
  v_ff.x <- xfs[t0+dt] - xfs[t0]
  v_ff.y <- yfs[t0+dt] - yfs[t0]
  v_ff.norm <- sqrt(v_ff.x^2 + v_ff.y^2)
  u_ff.x <- v_ff.x / v_ff.norm
  u_ff.y <- v_ff.y / v_ff.norm
  
  #get unit vector u_fl from follower current position to leader current position
  v_fl.x <- xs[leader, t0] - xfs[t0]
  v_fl.y <- ys[leader, t0] - yfs[t0]
  v_fl.norm <- sqrt(v_fl.x^2 + v_fl.y^2)
  u_fl.x <- v_fl.x / v_fl.norm
  u_fl.y <- v_fl.y / v_fl.norm
  
  #get projection of u_ff on vector perpendicular to group heading, proj_ff
  proj_ff <- u_ff.x * head.perp.x + u_ff.y * head.perp.y
  
  #get projection of u_fl on vector perpendicular to group heading, proj_fl
  proj_fl <- u_fl.x * head.perp.x + u_fl.y * head.perp.y
  
  #get overall turning leadership score
  score <- proj_ff * sign(proj_fl)
  
  #return turning leadership score
  return(score) 
}

#Get positions relative to centroid position and heading
#INPUTS:
# xs: [n.inds x n.times matrix] of eastings of each individual at each time point
# ys: [n.inds x n.times matrix] of northings of each individual at each time point
# centr.x: [n.times x 1 vector] of group centroid easting
# centr.y: [n.times x 1 vector] of group centroid northing
# headings: [n.times x 1 vector] of heading angles of the group (in radians)
#OUTPUTS:
# out: an object containing
# $xs.rel [n.inds x n.times matrix] of x position relative to centroid (positive x axis points in direction of group travel)
# $ys.rel [n.inds x n.times matrix] of y position relative to centroid (positive y axis points left relative to direction of travel)
get.rel.positions <- function(xs, ys, centr.x, centr.y, headings){
  
  #get displacement along eastings and northings of each individual from the centroid
  dxs <- xs - matrix(rep(centr.x, each = nrow(xs)), nrow = nrow(xs), ncol = ncol(xs))
  dys <- ys - matrix(rep(centr.y, each = nrow(ys)), nrow = nrow(xs), ncol = ncol(xs))
  
  #set up matrix of group headings
  angs <- matrix(rep(headings, each = nrow(ys)), nrow = nrow(xs), ncol = ncol(xs))
  
  #rotate displacements by -angs to get into group reference frame (direction of motion is along x axis)
  theta <- -angs
  xs.rel <- cos(theta)*dxs - sin(theta)*dys
  ys.rel <- sin(theta)*dxs + cos(theta)*dys
  
  out <- list()
  out$xs.rel <- xs.rel
  out$ys.rel <- ys.rel
  
  return(out)
  
}

#For a set of paired data (xvals, yvals), make a plot by binning data based on xvals and then plotting the mean (or another test statistic) of yvals within each bin
#INPUTS:
# xvals: vector of x values (used for binning)
# yvals: vector of y values (we will look at the test statistic of these within each bin)
# xbins: vector of x bins
# test.stat: which test statistic to use ('mean' and 'median' currently supported, defaults to 'mean')
# error.bar.quantiles: [string or numeric] if numeric (0-100), the associated interval of the distribution of yvalues for each bin will be plotted.
#     if 'boot', bootstrap confidence intervals will be plotted (1000 bootstraps, 95% CIs, assumes independence of points)
#     if 'clopper', the clopper pearson intervals will be plotted (only makes sense for binary data)
# xlab, ylab: [string] axis labels
# show.plot: [boolean] whether to show plot
# min.data.per.bin: [numeric] minimum number of data points needed to plot results (otherwise will be replaced with NA), defaults to 0
# plot.on.image: [boolean] if true, an alternative plot showing the probability distribution of y conditioned on x will be plotted (using a kernel density estimate)
#OUTPUTS:
# out: object of output containing
#   $xbins: x bins
#   $vals: test statistics of y values at each bin
#   $tots: total number of data points in each bin (non-NA)
bin.plot <- function(xvals, yvals, xbins, test.stat = 'mean', error.bar.quantiles = NULL, xlab='', ylab = '', show.plot = T, min.data.per.bin = 0, plot.on.image = F, h=NULL){
  
  test.stats <- uppers <- lowers <- tots <- rep(NA, length(xbins)-1)
  if(is.numeric(error.bar.quantiles)){
    pctile <- (100 - error.bar.quantiles)/2
  }
  for(i in 1:(length(xbins)-1)){
    vals.curr <- yvals[which((xvals >= xbins[i]) & (xvals < xbins[i+1]))]
    tots[i] <- sum(!is.na(vals.curr))
    if(length(vals.curr) > min.data.per.bin){
      if(test.stat=='mean'){
        test.stats[i] <- mean(vals.curr,na.rm=T)
      }
      if(test.stat == 'median'){
        test.stats[i] <- median(vals.curr, na.rm=T)
      }
      
      if(!is.null(error.bar.quantiles)){
        if(error.bar.quantiles == 'boot'){
          nboots <- 1000
          boots <- rep(NA, nboots)
          for(j in 1:nboots){
            if(test.stat=='mean'){
              boots[j] <- mean(sample(vals.curr, replace=T) , na.rm=T)
            }
            if(test.stat == 'median'){
              boots[j] <- median(sample(vals.curr, replace=T) , na.rm=T)
            }
          }
          lowers[i] <- quantile(boots, .025, na.rm=T)
          uppers[i] <- quantile(boots, .975, na.rm=T)
        } else if(error.bar.quantiles == 'clopper'){
          binom <- binom.test(sum(vals.curr==T,na.rm=T), sum(!is.na(vals.curr)))
          uppers[i] <- binom$conf.int[2]
          lowers[i] <- binom$conf.int[1]
        } else{
          lowers[i] <- quantile(vals.curr, pctile/100, na.rm=T)
          uppers[i] <- quantile(vals.curr, 1 - pctile/100, na.rm=T)
        }
      }
    }
    
  }
  
  out <- list()
  out$xbins <- xbins
  out$vals <- test.stats
  out$tots <- tots
  
  if(!is.null(error.bar.quantiles)){
    out$uppers <- uppers
    out$lowers <- lowers
  }
  
  if(show.plot & plot.on.image){
    mids <- xbins[1:(length(xbins)-1)]+diff(xbins)/2
    if(!is.null(h)){
      im <- kde2d(x = xvals, y = yvals, n = 100, lims = c(min(xbins,na.rm=T), max(xbins,na.rm=T), min(yvals), max(yvals)), h = h)
    } else{
      im <- kde2d(x = xvals, y = yvals, n = 100, lims = c(min(xbins,na.rm=T), max(xbins,na.rm=T), min(yvals), max(yvals)))
    }
    image.plot(im$z / rowSums(im$z), x = im$x, y = im$y, col = viridis(256), xlab = xlab, ylab = ylab)
    lines(mids, test.stats, lwd = 2, col = 'white')
    points(mids, test.stats, pch = 19, cex = 2, col = 'white')
  }
  
  if(show.plot & !plot.on.image){
    if(!is.null(error.bar.quantiles)){
      plot(NULL, xlim=c(min(xbins,na.rm=T),max(xbins,na.rm=T)), ylim = c(min(lowers,na.rm=T), max(uppers,na.rm=T)), xlab=xlab, ylab=ylab)
    } else{
      plot(NULL, xlim=c(min(xbins,na.rm=T),max(xbins,na.rm=T)), ylim = c(min(test.stats,na.rm=T), max(test.stats,na.rm=T)),xlab=xlab, ylab=ylab)
    }
    mids <- xbins[1:(length(xbins)-1)]+diff(xbins)/2
    if(!is.null(error.bar.quantiles)){
      polygon(x = c(mids, rev(mids)), y = c(uppers, rev(lowers)), border = NA, col = '#00000033')
    }
    lines(mids, test.stats)
    points(mids, test.stats, pch = 19)
  }
  
  return(out)
}






