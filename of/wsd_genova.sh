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

saveRDS(zz,'/home/matteo/own_data/PoD/topics/aedes_genmod/data/t_genova.RDS')