##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
## R code developed by: 
## Matteo Marcantonio
## ---------------------------------------------------##
## Latest update: September 2020
## ---------------------------------------------------##
## Scope: Visualise DynamAedes.R outputs.
## Must be run after run_DynamAedes.R.
## simout is the output of DynamAedes.R.
## ---------------------------------------------------##
## DOI: 
##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%##
### Load required packages
library(ggplot2)
library(raster)
library(sp)
library(gstat)
library(spatstat)
library(maptools)
library(rgeos)
library(geosphere)

#%%%%%%%%%%%%%%%%%%%#
### Define some variables
## Derive maxium simulation length in days from run_DynamAedes.R output 
days <- max(sapply(simout, function(x) length(x)))
## Define the day of introduction
str = 121
## Define the end-day of life cycle
endr = 365*2
## Define temperature dates, assuming the two years are 2019 and 2020
tdates <- seq(as.Date("2020-01-01"),as.Date("2020-01-01")+ndays,by="day")

#%%%%%%%%%%%%%%%%%%%##
### Derive probability of a successfull introduction at the end of the simulated period
pofe <- sum(unlist(lapply(simout, function(x) {
	days <- length(x)
	pe <- length(which(sum(x[days][[1]])>0)) 
	gc()
	return(pe)
}))) / length(simout)

### Print the probability
cat("\n### p of successfull introduction is:",pofe," ###\n")

#%%%%%%%%%%%%%%%%%%%#
## Derive abundance 95% CI for each life stage and in each day
# Define function
dabu95ci <- function(outl=NA,st=1,cores=1,days=0){
	out <- apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=c(0.025,0.50,0.975),na.rm=T);
	colnames(out)<-NULL
	outo <- rbind.data.frame(out,
		stage=rep(st,nrow(out)),
		day=as.factor(1:ncol(out)))
	return(t(outo))
}
# Apply function and format dataset for ggplot
dabu_df <- rbind.data.frame(dabu95ci(simout,1,days=days),dabu95ci(simout,2,days=days),dabu95ci(simout,3,days=days))
dabu_df$stage <- as.factor(dabu_df$stage)
levels(dabu_df$stage) <- c("Immature","Egg","Adult")
dabu_df$date <- rep(tdates[(str+1):(length(tdates))],3)
# Plot simple abundance time trend
ggplot(dabu_df, aes(y=`50%`+1,x=date,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`2.5%`+1,ymax=`97.5%`+1,fill=stage),col="white",alpha=0.1,outline.type="full") +
geom_smooth(linetype=1,size=1.5,se=FALSE) + labs(y="log10(Abundance)") +
scale_y_continuous(trans="log10") +
xlab("Date") +
theme(legend.pos="top") +
ggtitle("Smoothed 95% CI abundance per stage")

#%%%%%%%%%%%%%%%%%%%##
### Derive abundance median and inter-quartile per cell per day
## `st` defines the stage; 1=adult, 2=immature, 3=egg
st = 1 
cabu75ci <- function(outl=NA,st=st,days=0,ppp=c(0.25,0.50,0.75)){
	out <- apply(do.call(rbind.data.frame,lapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})})),2,quantile,probs=ppp,na.rm=T);
	colnames(out)<-NULL
	outo <- rbind.data.frame(out,
		stage <- rep(st,nrow(out)),
		day <- as.factor(1:ncol(out)))
	return(t(outo))
}
## Apply the function
cabu <- rbind.data.frame(cabu75ci(simout,st=1,days=days),cabu75ci(simout,st=2,days=days),cabu75ci(simout,st=3,days=days))
## Format columns
names(cabu) <- c("a25%","a50%","a75%","stage", "date")
cabu$stage <- as.factor(cabu$stage)
levels(cabu$stage) <- c("Egg","Immature","Adult")
cabu$date <- rep(tdates[(str+1):length(tdates)],3)

#%%%%%%%%%%%%%%%%%%%#
### Derive median and interquartile number of invaded cells per day
icel75ci <- function(outl=NA,days=0,ppp=c(0.25,0.50,0.75)){
	out <- apply(do.call(rbind.data.frame,lapply(outl, function(x) {
		lapply(1:days, function(y) {
			if( y<=length(x) ) {
				length(which(x[[y]]>0))
			} else {NA}})})),2,quantile,probs=ppp,na.rm=T)
	colnames(out)<- NULL
	outo <- rbind.data.frame(out,
		day=as.factor(1:ncol(out)))
	return(t(outo))
}
### Apply function
icel <- rbind.data.frame(icel75ci(simout,days=days))
icel$day <- tail(tdates,-str)
names(icel) <- c("ic25%","ic50%","ic75%","date")

#%%%%%%%%%%%%%%%%%%%#
### Merge female abundance per day and number of invaded cells and add simulated temperature
dfp <- merge(cabu[cabu$stage%in%"Adult",], icel, by="date")
dfp$temp <- colMeans(mat)[(str):endr]
dfp[,c(6:8)] <- dfp[,c(6:8)]*res/1000
### Plot the three trends together
ggplot(dfp, aes(y=`a50%`/`ic50%`, x=date, col="Adult females")) + 
geom_line(aes(y=temp*5, col="Mean daily temperature"),linetype=1) + 
geom_line(aes(y=`ic50%`*5, col="Invaded area (km^2)"),linetype=1) +
geom_ribbon(aes(ymin=`ic25%`*5,ymax=`ic75%`*5), col="transparent",alpha=0.4,show.legend=FALSE) +
geom_line(linetype=1,size=0.5,se=FALSE,alpha=1) +
scale_color_manual(values = c("#1B9E77","black","#A6CEE3","grey60")) +
scale_fill_manual(values = c("#1B9E77","black","#A6CEE3","grey60")) +
labs(x=NULL,y=expression(paste("Adult female abundance ", (km^-1)))) +
guides(colour=guide_legend(title="")) +
scale_y_continuous(sec.axis = sec_axis(name=expression(paste("Temperature (°C) -- Invaded area ", (km^2))),trans= ~ ./5), expand=c(0,0) ) +
scale_x_date(date_breaks = "3 months", date_labels = "%b-%y",expand=c(0.01,0.01)) +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16)) +
ggtitle("Median female abundance and interquartile Invaded area")

#%%%%%%%%%%%%%%%%%%%#
### Derive Eucliden distance from each cell introduction for all cells
## Define function
eudis <- function(x,x1,y,y1) {
	d <- sqrt((x1-x)^2 + (y1-y)^2)
	return(d)
}
## Select cells of introduction
introcells <- sapply(simout, function(x) which(colSums(x[[1]])>0))
### Apply the function to derive distances
coordd <- lapply(introcells, function(x) {
	ccc <- apply(cc, 1, function(y) {
		ccc <- eudis(x=cc[unlist(x),1],x1=y[1],y=cc[unlist(x),2],y1=y[2])
	})
	cc$distm <- ccc
	return(cc)
})

#%%%%%%%%%%%%%%%%%%%#
### Select only cells with at least 1 mosquito and match cells with distance from cell from introduction
ainvc <- lapply(1:length(simout), function(d) {
	sapply(simout[[d]], function(y) {
		z <- data.frame(y[,which(colSums(y)>0)])
		colnames(z) <- which(colSums(y)>0)
		z[4,] <- coordd[[d]]$distm[which(colSums(y)>0)]
		return(z)
	})
})

#%%%%%%%%%%%%%%%%%%%#
### Derive interquartile of distances from cell of introduction for each day
interdis <- function(outl=NA,days=0){
	out <- do.call(rbind.data.frame, lapply(1:days, 
		function(x) {
			quantile(unlist(sapply(1:length(outl),
				function(y) {
					mean(if(x<length(outl[[y]])) as.integer(outl[[y]][[x]][4,]) else NA,na.rm=T)
				})), probs=c(0.25,0.50,0.75), na.rm=T)
		}))
	names(out) <- c("25%","50%","75%")
	return(out)
}
### Apply the function
## Format date
asp <- interdis(ainvc,days=endr)
asp$date <- head(tdates,-1)
### Plot the interquartile of distances at which mosquitoes are found every day
ggplot(asp, aes(y=`50%`,x=date)) + 
geom_line() +
geom_ribbon(aes(ymin=`25%`,ymax=`75%`),colour="grey90",alpha=0.4) +
geom_smooth(se=FALSE) +
scale_color_brewer(palette = "Dark2") +
scale_fill_brewer(palette = "Dark2") +
labs(x=NULL,y="Distance (km) from cell of introduction") +
scale_y_continuous(labels=function(x) format(x/1000, scientific = FALSE)) +
scale_x_date(date_breaks = "3 months", date_labels = "%b\n%y",expand=c(0.01,0.01)) +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))+
ggtitle("Spread of the invasive mosquito population")

#%%%%%%%%%%%%%%%%%%%#
### Some bonus plots and information
## Plot crude trend in space for one single iteration
# Select simulation to visualise
simu <- 2
# Reduce over time
summ <- Reduce(`+`,simout[[simu]])
# Use `r` raster defined in run_DynamAedes.R as template 
rsum <- r
# Plug in crude average abundace per day per cell
values(rsum)<- summ[3,]/days
# Define an adequate color palette
maxColorValue <- max(values(rsum))
palette <- colorRampPalette(c("white","blue","red"),bias=100)(maxColorValue)
# Plot invasion trend in space with road segment
raster::plot(rsum,col = palette)
plot(buff,add=T,col = rgb(red = 1, green = 1, blue = 0, alpha = 0.5))
title("Crude average abundance of female mosquitoes\nalong the simulation period\n& road segment")

## Get date and median number of adults at population peak
dfp[dfp$stage%in%"Adult",][which.max(dfp[dfp$stage%in%"Adult",]$`a50%`),]

## Get date and maximum interquartile distance travelled by mosquitoes
asp[which.max(asp$`75%`),]