#get dailty temperature estimates
require(raster)
require(microclima)
require(NicheMapR)

setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/venezia/')
# venezia<-raster('venezia_region3035.tiff' )
dem250<-get_dem(venezia,resolution=250,zmin=0)
# rplot(venezia)
venezia250<-raster("venezia_dtm_clipped.tif")

#Add precipitation
prcp <- wda[[1]]$prcp[wda[[1]]$yearmoda>="2017-01-01"&wda[[1]]$yearmoda<"2019-01-01"]
prcp <- ifelse(is.na(prcp),0,prcp)

ventemps <- microclima::runauto(venezia250, dstart="01/01/2017", dfinish="31/12/2018", hgt = 2,	l = NA, x = NA, r.is.dem=TRUE, habitat = 14, plot.progress = FALSE, save.memory=TRUE, coastal=FALSE, use.raster=TRUE, zmin=0, summarydata = FALSE, dailyprecip=prcp)

#get daily mean
day<-seq(1, dim(ventemps$temps)[3], 24)
dailyT<-list()

# outr <- lapply(ventemps$temps, function(x) {
# 	mclapply(seq(1, dim(ventemps$temps)[3], 24), function(i) {
# 		cat("day: ",(i+(23/24)))
# 		tday <- apply(x[,,i:(i+23)], c(1,2), FUN = mean)
# 		#tdayr <- raster(tday,template=dem1kmr)
# 		return(tday)
# 	}, mc.cores=47)
# })

for(i in day){
	j <- i + 23
	print(j/24)
	Tday <- apply(ventemps$temps[,,i:j], c(1,2), FUN = mean )
	dailyT <- append(dailyT,raster(Tday,template=venezia250))
}

lname<-as.character(paste("d", "_", seq(1:(dim(ventemps$temps)[3]/24)), sep=""))
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
saveRDS(T_df, '~/venezia_data250mno_coast_prcp_2016.RDS')
#write.csv(T_df, '/home/matteo/GitHub/euaeae/t_data.csv')          