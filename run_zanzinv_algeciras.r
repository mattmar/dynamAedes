## Load time series of daily average temperature and distance matrices
setwd("~/GitHub/euaeae")
#Start from 1st of January 2017
w <- readRDS('data/algeciras_temps11-19_corr.RDS')[,c(1:3,2193:3229)] #select years 1=2011; 365*6 (2193) means you start from 1 Jan 2017
ww <- w[,-c(1)] #remove index
cc <- ww[,c(1:2)] #coordinates
ww <- ww[,-c(1:2)] #remove lat long
pld <- read.csv('data/algeciras_distmatrix.txt')
row.names(pld) <- which(w[,1]%in%pld[,1])
colnames(pld) <- row.names(pld)
pld <- apply(pld,2,function(x) round(x/10000,1)*10000) #make round 100's 
gc()

#intro_cells for algeciras; this is the index and not the rowname
introvector<-c(15664,15665,15666,15667,15685,15686,15687,15688,15704,15705,15706,15707,15708,15723,15724,15725,15726,15727,15741,15742,15743,15744,15745,15757,15758,15759,15760,15761,15762,15774,15775,15776,15777,15778,15779,15780,15792,15793,15794,15795,15796,15797,15798,15810,15811,15812,15813,15814,15815,15816,15832,15833,15834,15835,15836,15837,15838,15854,15855,15856,15857,15858,15859,15860,15875,15876,15877,15878,15879,15880,15881,15894,15895,15896,15897,15898,15899,15900,15912,15913,15914,15915,15916,15917,15929,15930,15933,15934)

# Run simulations, start 15th of July 2017; end 31 October 2019
str=235
endr=1037
source("zanzinv.r")
zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld,startd=str,endd=endr,n.clusters=47, cluster.type="SOCK",iter=100,intro.cells=introvector,intro.eggs=250,sparse.output=FALSE,compressed.output=TRUE)

str=225
endr=445
zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld,startd=str,endd=endr,n.clusters=7, cluster.type="SOCK",iter=7,intro.cells=introvector,intro.eggs=100,sparse.output=FALSE,compressed.output=TRUE)

#temps.matrix=ww; cells.coords=cc; road.dist.matrix=pld;startd=210; endd=240; n.clusters=2; cluster.type="SOCK" ; iter=2; intro.cell=NA; intro.adults=100

### Plot results ###
#trend in time, just one iteration
outlp<-mclapply(zanzout[[5]], function(x) {if(!is.na(x[[1]])) apply(x, MARGIN=c(1, 2), sum) else NULL},mc.cores=4)
outlp<-outlp[!sapply(outlp, function(x) is.null(x[[1]]))]
outlt<-melt(t(sapply(outlp, function(x) rowSums(x))))
outlt$day<-as.integer(row.names(outlt))
ggplot(outlt, aes(x=Var1,y=value,col=as.factor(Var2))) + geom_line()

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
citoplot$date<-as.Date(citoplot$day,origin="2017-07-25")

w <- readRDS('~/GitHub/euaeae/data/algeciras_temps11-19_corr.RDS')[,2193:3229][,str:((str-1)+days)]

citoplot$tempm<-rep(as.numeric(w[15635,]),3)

citoplot$tempm[which(citoplot$tempm>40000)]<-NA

g1<-ggplot(citoplot, aes(y=`50%`,x=date,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=stage),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Populazion size") +
geom_line(aes(y=tempm/10),col="black") + labs(y="Populazion size") +
ggtitle("Introduction of 100 eggs in Algeciras harbour") +
scale_y_continuous(sec.axis = ~ ./100) +
xlim(range(seq(as.Date("2017-01-01")+str, by = "day", length.out = 100))); g1

ggsave("~/algeciras_500.png",g1,dpi=400,scale=1.5,width=25,height=15,unit="cm")

###Spatial spread
w <- readRDS('data/algeciras_temps11-19_corr.RDS')[,c(1:3,2193:3229)] #select years 1=2011;
cc <- w[,c(2:3)] #coordinates
ccc <- coordinates(spTransform(SpatialPoints(coords = cc, proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))
library(geosphere)
cc$distm<-as.integer(distm(ccc[15664,],ccc))

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

plotdistt<-cbind.data.frame(returndis(zanzred,days=days),day=1:days)
plotdistt$date<-as.Date(plotdistt$day,origin="2017-07-25")

g2 <- ggplot(plotdistt, aes(y=`50%`,x=date)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Distance (m)") +
ggtitle("Introduction of 100 eggs in Venice harbour"); g2

ggsave("~/distance_travelled_algeciras_500.png",g2,dpi=400,scale=1.5,width=25,height=15,unit="cm")

#trend in space, just one iteration
library(rgdal)
wadults<-cbind(w[2:3],sapply(outlp, function(x) x[3,]))
wlarvae<-cbind(w[2:3],sapply(outlp, function(x) x[2,]))
weggs<-cbind(w[2:3],sapply(outlp, function(x) x[1,]))
roads<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/algeciras/algeciras_roads.sqlite')
area<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/algeciras/algeciras_urban.sqlite')

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
