#get dailty temperature estimates
require(raster)
require(microclima)
require(NicheMapR)

setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/genova/')
dem1km<-raster('dem1k_clipped.tiff')
dem1kmr<-projectRaster(dem1km, crs = "+init=epsg:3395 +proj=merc +lon_0=0 +k=1 +x_0=0 +y_0=0 +datum=WGS84" )
nedem100<-get_dem(dem1kmr,resolution=100,zmin=1)
nedem250<-get_dem(dem1kmr,resolution=250,zmin=1)
nedem1km<-get_dem(dem1kmr,resolution=1000,zmin=1)


#Define period
#tme <- as.POSIXlt(c(0:(24*365)) * 3600, origin ="2016-01-01 00:00", tz = "UTC")
#Donwload NCE
#ho<-hourlyNCEP(lat=44.407089, long= 8.954859, tme=tme, reanalysis2 = TRUE)
#coste<-coastalNCEP(dem1kmr, ho, steps = 8, use.raster = T, zmin = 0,plot.progress = T, tidyr = F)
#md<-micro_ncep(loc = c(44.407089, 8.954859), dstart = "15/01/2015", dfinish = "23/01/2015", REFL = 0.15, slope = 0, aspect = 0, DEP = c(0, 2.5,  5,  10,  15,  20,  30,  50,  100,  200), minshade = 0, maxshade = 90, Usrhyt = 0.01, coastal=TRUE,hourlydata=ho)
prcp <- data$prcp[data$yearmoda>="2016-01-01"&data$yearmoda<"2017-01-01"]

temps <- microclima::runauto(nedem100, dstart="01/01/2016", dfinish="31/12/2016", hgt = 1,	l = NA, x = NA, r.is.dem=TRUE, habitat = 14, plot.progress = FALSE, save.memory=TRUE, dailyprecip=prcp,coastal=TRUE,use.raster=TRUE,zmin=1)

#get daily mean
day<-seq(1, dim(temps$temps)[3], 24 )
dailyT<-list()

outr <- lapply(temps$temps, function(x) {
	mclapply(seq(1, dim(temps$temps)[3], 24), function(i) {
		cat("day: ",(i+(23/24)))
		tday <- apply(x[,,i:(i+23)], c(1,2), FUN = mean)
		tdayr <- raster(tday,template=dem1kmr)
		return(tdayr)
	}, mc.cores=47)
})

for(i in day){
	j <- i + 23
	print(j/24)
	Tday <- apply(temps$temps[,,i:j], c(1,2), FUN = mean )
	dailyT <- append(dailyT,resample(raster(Tday,template=nedem),nedem250, method='bilinear'))
}

lname<-as.character(paste("d", "_", seq(1:(dim(temps$temps)[3]/24)), sep=""))
names(dailyT)<-lname

length(dailyT)
str(dailyT)
class(dailyT)
dailyT

rtest<-dailyT$d_5
raster::plot(rtest)
# saveRDS(dailyT, "dailyT.rds")

#extract temperature TS for each pixel
pts<-data.frame(coordinates(raster(dailyT$d_1)))
coordinates(pts)<-~x+y
raster::plot(rtest)
raster::plot(pts, add=T)

T_df<-lapply(dailyT, function(x) {y=extract(x, pts);return(y)})
T_df<-as.data.frame(c(pts,T_df))
summary(T_df)
T_df<-na.omit(T_df)
T_df[,1:2] <- coordinates(spTransform(SpatialPoints(coords = T_df[,c(1,2)], proj4string = CRS(projection(dailyT$d_1))), CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))

#T_df<-apply(T_df,2,function(x) as.integer(x*10))
write.csv(T_df, '~/t_data250m.csv')
#write.csv(T_df, '/home/matteo/GitHub/euaeae/t_data.csv')          