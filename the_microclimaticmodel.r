#get dailty temperature estimates
require(raster)
require(microclima)
require(NicheMapR)

setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/')
# venezia<-raster('venezia_region3035.tiff' )
# alg250<-get_dem(alg,resolution=250,zmin=0)
# rplot(venezia)
alg250<-raster("algeciras/alg_dtm_clipped.tif")
ven250<-raster("venezia/venezia_dtm.tif")

#Dates
sdate<-as.Date("2019-01-01")
edate<-as.Date("2019-10-31")

#Add precipitation
# prcp <- wda[[1]][wda[[1]]$yearmoda>=sdate&wda[[1]]$yearmoda<edate,c("yearmoda","prcp")]
# dd <- data.frame(yearmoda=seq(as.Date(sdate),as.Date(edate),1))
# prcp <- merge(dd,prcp,by="yearmoda",all.x=T)
# prcp <- ifelse(is.na(prcp$prcp),0,prcp$prcp)

ventemps <- microclima::runauto(ven250, dstart=format(sdate,"%d/%m/%Y"), dfinish=format(edate,"%d/%m/%Y"), hgt = 2, l = NA, x = NA, r.is.dem=TRUE, habitat = 14, plot.progress = FALSE, save.memory=TRUE, coastal=FALSE, use.raster=TRUE, zmin=0, summarydata = FALSE)

#get daily mean
day<-seq(1, dim(algtemps$temps)[3], 24)
dailyT<-list()

for(i in day){
	j <- i + 23
	print(j/24)
	Tday <- apply(algtemps$temps[,,i:j], c(1,2), FUN = mean )
	dailyT <- append(dailyT,raster(Tday,template=alg250))
}

lname<-as.character(paste("d", "_", seq(1:(dim(algtemps$temps)[3]/24)), sep=""))
names(dailyT)<-lname

#extract temperature TS for each pixel
pts<-data.frame(coordinates(raster(dailyT$d_1)))
coordinates(pts)<-~x+y

T_df<-lapply(dailyT, function(x) {y=extract(x, pts);return(y)})
T_df<-as.data.frame(c(pts,T_df))
summary(T_df)
T_df<-na.omit(T_df)
T_df[,1:2] <- coordinates(spTransform(SpatialPoints(coords = T_df[,c(1,2)], proj4string = CRS(projection(dailyT$d_1))), CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")))

saveRDS(T_df, '~/algeciras_data250mcoast_noprcp_2017.RDS')


# data<-as.data.frame(matrix(c(395,157,13,27,163,0,2,0,147,0,4,0,14,35,55,170,10,32,23,185,10,28,2,36,9,20,4,32,7,0,10,0,6,0,4,0,0,15,0,53,0,11,0,23,0,8,0,33,0,6,0,14),ncol=4,nrow=13,byrow=T))
# names(data) <-c("Ae.ae_posres","Ae.al_posres","Ae.ae_pospub","Ae.al_pospub")
# data$habitat <-c("Domestic containers", "Ornamental
# Containers", "Flower pot plate/tray","Drains", "Discarded
# receptacles","Plants","CanvasPlasticSheet","Puddle/ Ground depression","Roof top/RoofGutters","Gully trap","Inspection chamber","Bins","Puddle/Ground depression")

# datam<-melt(data,id.vars="habitat")
# datam$species<-sapply(strsplit(as.character(datam$variable),"_",fixed=FALSE), function(x) x[1])
# datam$loc<-sapply(strsplit(as.character(datam$variable),"_",fixed=FALSE), function(x) x[2])

# ggplot(data=datam, aes(x=habitat, y=value,fill=habitat,group=species)) + 
#   geom_bar(stat="identity") +
#   facet_wrap(~loc+species)