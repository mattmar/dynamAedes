#get dailty temperature estimates
require(raster)
require(microclima)
require(NicheMapR)

#lets do a test for genova
setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/genova/')
ge.shp<-shapefile("genova.shp")
ext<-extent(ge.shp)
# resol<-c(0.0008,0.0008)

aoi = raster(ext)
crs(aoi)<-"+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 "
res(aoi)=c(0.0083,0.0083)
aoi[]<-1
aoi
# aoi<-raster(res=resol, ext=ext, vals=1)

# dem500m <- get_dem(r=aoi, resolution = 500)
# dem200 <- get_dem(lat = 8.94, long = 44.40, resolution = 200)
dem1k <- get_dem(aoi, resolution = 0.0083)
writeRaster(dem1k, "dem1k", overwrite=T)
# plot(dem50m, add=T)
# plot(aoi)
# dem200<-projectRaster(dem200, crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 " )


# dem5m <- get_dem(lat = 49.97, long = -5.22, resolution = 5)
# r <- microclima::get_dem(lat = 38.467429, long = -28.398995, resolution = 30)


# Takes ~ c. 5 minutes to run
temps <- microclima::runauto(dem1k, "1/06/2018", "30/06/2018", hgt = 0.1,
                             l = NA, x = NA, coastal = F,
                             habitat = "Barren or sparsely vegetated",
                             plot.progress = TRUE)


str(temps$temps)
class(temps$temps)


#get daily mean
day<-seq(1, 720, 24 )
j<-vector()
# dailyT<-array(c(NA,NA), dim=c(17,52,1))
dailyT<-list()

for(i in day){
  a<-i+23
  j<-c(j, a)
  #print(j)
  Tday<-apply(temps$temps[,,i:j], c(1,2), FUN = mean )
  dailyT<- append(dailyT,raster(Tday))
}

lname<-as.character(paste("d", "_", seq(1:30), sep=""))
names(dailyT)<-lname

length(dailyT)
str(dailyT)
class(dailyT)
dailyT

rtest<-dailyT$d_5
plot(rtest)
# saveRDS(dailyT, "dailyT.rds")

#extract temperature TS for each pixel
pts<-data.frame(coordinates(raster(dailyT$d_1)))
coordinates(pts)<-~x+y
plot(rtest)
plot(pts, add=T)

T_df<-lapply(dailyT, function(x) {y=extract(x, pts);return(y)})
T_df<-as.data.frame(c(pts,T_df))
summary(T_df)
T_df<-na.omit(T_df)
# write.csv(T_df, "T_df.csv")
                    