###Download data from NOAA
##Genova Sestri is the weather station of reference whereas Albengs is just to use as a comparison
##This bits of code must be run in a Unix Shell (Bash)
cd '~/own_data/PoD/topics/aegypti_eu/spatial_data/genova/temperature/'
for year in `seq 2010 2019`; do wget -P gs_ws ftp://ftp.ncdc.noaa.gov/pub/data/gsod/$year/161200-99999-$year.op.gz ;done

##R now
#load data
setwd('~/own_data/PoD/topics/aegypti_eu/spatial_data/genova/temperature/gs_ws/')
dgs <- do.call(rbind.data.frame,lapply(list.files(pattern="*.gz")[grep(".gz",list.files(pattern="*.gz"))], read.table,sep="",skip=1,header=FALSE))[,c(1,3,4,20)]
names(dgs) <- c("id","yearmoda","temp","prcp")

#Transform in a list and process it
wda <- list(dgs)

#Fix units and remove NAs using NOAA labels
wda <- lapply(wda, function(x) {
	x$yearmoda <- as.Date(as.character(x$yearmoda),"%Y%m%d"); 
	x$temp <- round((x$temp-32)/9*5,1); 
	x$prcp[which(x$prcp%in%"99.99")] <- 0;
	x$prcp <- round(as.numeric(gsub("[A-Z]","",x$prcp))*25.4,2); return(x)
})

#Basic plots
plot(wda[[1]]$yearmoda,wda[[1]]$temp)
plot(wda[[1]]$yearmoda,wda[[1]]$prcp)

##Checks
#Format date and check monthly means
wda[[1]]$month<-format(wda[[1]]$yearmoda,"%m")
aggregate(temp~month,wda[[1]],"mean")

#Check yearly means
wda[[1]]$year <- format(wda[[1]]$yearmoda,"%y")
mean(aggregate(prcp~year,wda[[1]],"sum")$prcp[3:10])
mean(aggregate(temp~year,wda[[1]],"mean")$temp)

#Interpolate lineraly to fill gaps in time series
library(zoo)
z <- read.zoo(wda[[1]][,c(2:3)])
zz <- data.frame(temp=na.approx(z, xout = seq(start(z), end(z), "day")))
zz$yearmoda <- as.Date(row.names(zz))

gentemp<-zz[zz$yearmoda>"2018-01-01",]$temp

#Load Microclima data
gen250 <- readRDS('~/own_data/PoD/topics/aegypti_eu/spatial_data/genova/temperature/genova_temps17-19.RDS')

#Extract row from raster 46666 for genova sestri weather station
tomerge <- cbind.data.frame(yearmoda=seq(as.Date("2017-01-01"),as.Date("2019-10-31"),1),temp=as.numeric(gen250[46666,-c(1:3)])/1000)

#Merge the two datasets
comp <- merge(tomerge, zz[,c("yearmoda","temp")],by="yearmoda",all.x=T)

#Use this vector of differences to correct the Microclima temperature just in the weather station pixel
comp$cgs<-comp$temp.y-comp$temp.x

#Plot the three time series (weather station, Microclima and difference) to visualize the mismatch
comp2 <- melt(comp,id.vars="yearmoda")
comp2$variable <- factor(comp2$variable)
levels(comp2$variable) <- c("Microclim250","weather station","diff")
ggplot(comp2,aes(x=yearmoda,y=value,color=variable,group=variable)) +geom_line() + labs(y="temperature C") +ggtitle("genova")

#Apply the correction to the full Microclima dataset (this is the input data then used in the model)
gen250c <- gen250
gen250c[,-c(1:3)] <- t(apply(gen250[,-c(1:3)],1, function(x) x+round(comp$cgs*1000,0)))
#saveRDS(gen250,'~/own_data/PoD/topics/aegypti_eu/spatial_data/genova/temperature/genova_temps17-19_corr.RDS')

## Visally compare original and corrected Microclima datasets with Albenga weather station data to see if the correction can be generalised to any pixel in the study area
#Load albenga weather station data
setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/genova/temperature/albenga/')
dgs <- do.call(rbind.data.frame,lapply(list.files(pattern="*.gz")[grep(".gz",list.files(pattern="*.gz"))], read.table,sep="",skip=1,header=FALSE))[,c(1,3,4,20)]

#Process it
names(dgs) <- c("id","yearmoda","temp","prcp")
wda <- list(dgs)

wda <- lapply(wda, function(x) {
	x$yearmoda <- as.Date(as.character(x$yearmoda),"%Y%m%d"); 
	x$temp <- round((x$temp-32)/9*5,1); 
	x$prcp[which(x$prcp%in%"99.99")] <- 0;
	x$prcp <- round(as.numeric(gsub("[A-Z]","",x$prcp))*25.4,2); return(x)
})

#Extract the row matching Albenga weather station geolocation (54507) in the Microclima corrected dataset
tomerge0 <- cbind.data.frame(yearmoda=seq(as.Date("2017-01-01"),as.Date("2019-10-31"),1),temp=as.numeric(gen250[54507,-c(1:3)])/1000)
tomerge1 <- cbind.data.frame(yearmoda=seq(as.Date("2017-01-01"),as.Date("2019-10-31"),1),temp=as.numeric(gen250c[54507,-c(1:3)])/1000)

comp <- merge(merge(tomerge0, wda[[1]][,c("yearmoda","temp")],by="yearmoda",all.x=T),tomerge1,by="yearmoda",all.x=T) # some weather station T are NA's

#Derive vector of differences for plotting purposes
comp$cgs <- comp$temp.y-comp$temp.x

comp2 <- melt(comp,id.vars="yearmoda")
comp2$variable <- factor(comp2$variable)
levels(comp2$variable) <- c("Microclim250","weather station","Microclim250_corrected","diff")
ggplot(comp2,aes(x=yearmoda,y=value,color=variable,group=variable)) +geom_line() + labs(y="temperature C") +ggtitle("Albenga") +xlim(as.Date("2017-01-01"),as.Date("2019-10-31")) + theme(legend.position="bottom")