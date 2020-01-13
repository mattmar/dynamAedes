## Load time series of daily average temperature and distance matrices
setwd("~/GitHub/euaeae")
w <- readRDS('data/venezia_temps11-18.RDS')[,c(1:3,2193:2925)] #select years 1=2011;
ww <- w[,-c(1)] #remove index
cc <- ww[,c(1:2)] #coordinates
ww <- ww[,-c(1:2)] #remove lat long
pld <- read.csv('data/venezia_distmatrix.txt')
row.names(pld) <- which(w[,1]%in%pld[,1])
colnames(pld) <- row.names(pld)
pld <- apply(pld,2,function(x) round(x/10000,1)*10000) #make round 100's 
gc()

#intro_cells for venice; this is the index and not the rowname
introvector<-c(27475,27476,27645,27646,27647,27648,27820,27821,27822,27823,27824,27992,27993,27994,27995,27996,28160,28161,28162,28163,28164,28323,28324,28325,28326,28327,28328,28329,28481,28482,28483,28484,28485,28486,28487,28488,28489,28490,28625,28626,28627,28628,28629,28630,28631,28632,28633,28634,28768,28769,28770,28771,28772,28773,28774,28775,28776,28777,28913,28914,28915,28916,28917,28918,28919,28920,28921,28922,28923,28924,29070,29071,29072,29073,29074,29075,29076,29077,29078,29079,29080,29081,29082,29083,29235,29236,29237,29238,29239,29240,29241,29242,29243,29244,29245,29246,29247,29248,29249,29250,29408,29409,29410,29411,29412,29413,29414,29415,29416,29417,29418,29419,29420,29584,29585,29586,29587,29588,29589,29590,29591,29592,29593,29594,29595,29596,29771,29772,29773,29774,29775,29776,29777,29778,29779,29780,29781,29782,29783,29784,29785,29980,29981,29982,29983,29984,29985,29986,29987,29988,29989,30197,30198,30199,30200,30201,30202,30203,30204,30205,30206,30207,30208,30418,30419,30420,30421,30422,30423,30424,30425,30426,30427,30428,30429,30430,30431,30641,30642,30643,30644,30645,30646,30647,30648,30649,30650,30651,30652,30653,30865,30866,30867,30868,30869,30870,30871,30872,30873,30874,30875,30876,30877,31097,31098,31099,31100,31101,31102,31103,31104,31105,31332,31333,31334,31335,31336,31337,31338,31339,31567,31568,31569,31570,31571,31572,31573,31574,31575,31810,31811,31812,31813,31814,31815,31816,31817,31818,31819,31820,31821,32063,32064,32065,32066,32067,32068,32069,32070,32071,32072,32073,32074,32075,32318,32319,32320,32321,32322,32323,32324,32325,32326,32327,32328,32577,32578,32579,32580,32581,32582,32583,32584,32585,32840,32841,32842,32843,32844,32845,32846,33096,33097,33098,33099,33100,33101,33353,33354)

# Run simulations
source("zanzinv.r")
zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld,startd=196,endd=733,n.clusters=47, cluster.type="SOCK",iter=100,intro.cells=introvector,intro.eggs=250,sparse.output=FALSE,compressed.output=TRUE)

#temps.matrix=ww; cells.coords=cc; road.dist.matrix=pld;startd=210; endd=240; n.clusters=2; cluster.type="SOCK" ; iter=2; intro.cell=NA; intro.adults=100

### Plot results ###
#trend in time, just one iteration
outlp<-mclapply(zanzout[[5]], function(x) {if(!is.na(x[[1]])) apply(x, MARGIN=c(1, 2), sum) else NULL},mc.cores=4)
outlp<-outlp[!sapply(outlp, function(x) is.null(x[[1]]))]
outlt<-melt(t(sapply(outlp, function(x) rowSums(x))))
outlt$day<-as.integer(row.names(outlt))
ggplot(outlt, aes(x=Var1,y=value,col=as.factor(Var2))) + geom_line()

#trend in space, just one iteration
library(rgdal)
wadults<-cbind(w[2:3],sapply(outlp, function(x) x[3,]))
wlarvae<-cbind(w[2:3],sapply(outlp, function(x) x[2,]))
weggs<-cbind(w[2:3],sapply(outlp, function(x) x[1,]))
roads<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/venezia/venezia_roads.sqlite')
area<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/venezia/venezia_urban.sqlite')

for (i in 1:(ncol(wadults)-2)) {
	png(paste('/home/matteo/own_data/PoD/topics/aegypti_eu/figures/figure_',sprintf("%03d", i),".png",sep=""),width = 480*3, height = 480*2)
	raster::plot(wadults[,1],wadults[,2],cex=0.1,col="black",type="p",xlab="",ylab="")
	raster::plot(wadults[,1],wadults[,2],cex=0.1,col="black",type="p",xlab="",ylab="")
	raster::plot(wadults[,1],wadults[,2],cex=0.1,col="black",type="p",xlab="",ylab="")
	raster::plot(area,add=T)
	raster::plot(roads,add=T)
	text(x=4212000,y=2380000,as.Date(i+150, origin = "2016-01-01"),col="red",cex=2)
	points(weggs[,1][which(weggs[,i+2]>0)],weggs[,2][which(weggs[,i+2]>0)],pch=16,cex=1.1,col="blue")
	points(wlarvae[,1][which(wlarvae[,i+2]>0)],wlarvae[,2][which(wlarvae[,i+2]>0)],pch=16,cex=1.1,col="green")
	points(wadults[,1][which(wadults[,i+2]>0)],wadults[,2][which(wadults[,i+2]>0)],pch=16,cex=1.1,col="red")
	#text(wadults[,1],wadults[,2],labels=ifelse(wadults[,i+2]>0,wadults[,i+2],""),cex=wadults[,i+2]/10,col="black",pos=2,offset=runif(1,0:1))
	text(x=4270000,y=2380000,labels=paste("Total number of E",sum(weggs[,i+2]),"L",sum(wlarvae[,i+2]),"A",sum(wadults[,i+2]),sep=": "),col="black",cex=2)
	#Sys.sleep(0.1)
	dev.off()
}	

#Derive number of time an iteration is still active at the end of next spring
library(parallel)
days<-max(sapply(zanzout,length))
#outlp<-mclapply(zanzout, function (x) {lapply(x, function(y) {if(!is.na(y[[1]])) apply(y, MARGIN=c(1, 2), sum) else NULL})},mc.cores=20)

#check number of individuals in each stage
#lapply(outlp, function(x) if(!is.null(x[days])) rowSums(x[days][[1]]))

#derive 95% CI in each day
returnci<-function(outl=NA,st=1,cores=1,days=0){
	out<-apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=c(0.025,0.50,0.975),na.rm=T);
	colnames(out)<-NULL
	outo<-rbind.data.frame(out,
		stage=rep(st,nrow(out)),
		day=as.factor(1:ncol(out)))
	return(t(outo))
}

citoplot<-rbind.data.frame(returnci(zanzout,1,days=days),returnci(zanzout,2,days=days),returnci(zanzout,3,days=days))

citoplot$stage<-as.factor(citoplot$stage)
citoplot$date<-as.Date(citoplot$day,origin="2017-07-15")

w <- readRDS('~/GitHub/euaeae/data/venezia_temps11-18_corr.RDS')[,c(1:3,2193:2925)][,196:(195+days)]

citoplot$tempm<-rep(as.numeric(w[26102,]),3)

g1<-ggplot(citoplot, aes(y=`50%`,x=date,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=stage),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Populazion size") +
geom_line(aes(y=tempm/10),col="black") + labs(y="Populazion size") +
ggtitle("Introduction of 100 eggs in Venezia harbour") +
scale_y_continuous(sec.axis = ~ ./100) +
xlim(as.Date("2018-01-01"),as.Date("2019-04-30")); g1

ggsave("~/venezia_500.png",g1,dpi=400,scale=1.5,width=25,height=15,unit="cm")

###Spatial spread
w <- readRDS('data/venezia_temps11-18.RDS')[,c(1:3,((365*5)+3):((365*8)+3))] #select years 1=2011;
cc <- w[,c(2:3)] #coordinates
ccc <- coordinates(spTransform(SpatialPoints(coords = cc, proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))
library(geosphere)
cc$distm<-as.integer(distm(ccc[30426,],ccc))

#save only cells with at least 1 guy
zanzred <- mclapply(zanzout, function(x) {sapply(x, function(y) {
	z <- data.frame(y[,which(colSums(y)>0)]); 
	colnames(z) <- which(colSums(y)>0)
	z[4,]<-cc$distm[which(colSums(y)>0)]; 
	return(z)
})},mc.cores=1)

#derive 95% CI in each day
returndis <- function(outl=NA,days=0){
	out <- do.call(rbind.data.frame, lapply(1:days, 
		function(x) {
			quantile(unlist(sapply(1:length(outl),
				function(y) {mean(if(x<length(outl[[y]])) as.integer(outl[[y]][[x]][4,]) else NA,na.rm=T)})), probs=c(0.025,0.50,0.975), na.rm=T)
		}))
	names(out) <- c("2.5%","50%","97.5%")
	return(out)
}

plotdistt<-cbind.data.frame(returndis(zanzred,days=200),day=1:200)
plotdistt$date<-as.Date(plotdistt$day,origin="2016-06-15")

g2 <- ggplot(plotdistt, aes(y=`50%`,x=date)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Distance (m)") +
ggtitle("Introduction of 100 eggs in Venice harbour"); g2
ggsave("~/distance_travelled_venezia_500.png",g2,dpi=400,scale=1.5,width=25,height=15,unit="cm")
