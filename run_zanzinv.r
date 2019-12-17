## Load time series of daily average temperature and distance matrices
setwd("~/GitHub/euaeae")
w <- readRDS('data/venezia_temps11-18.RDS')
ww <- w[,-c(1)] #remove index
cc <- ww[,c(1:2)] #coordinates
ww <- ww[,-c(1:2)] #remove lat long
pld <- read.csv('data/venezia_distmatrix.txt')
row.names(pld) <- which(w[,1]%in%pld[,1])
colnames(pld) <- row.names(pld)
pld <- apply(pld,2,function(x) round(x/10000,1)*10000) #make round 100's 
gc()

# Run simulations
source("zanzinv.r")
zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld,startd=120,endd=515,n.clusters=40, cluster.type="SOCK",iter=1000,intro.cell=28486,intro.eggs=100,sparse.output=FALSE,compressed.output=TRUE)

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
lapply(outlp, function(x) if(!is.null(x[days])) rowSums(x[days][[1]]))

#derive 95% CI in each day
returnci<-function(outl=NA,st=1,cores=8,days=0){
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

g1<-ggplot(citoplot, aes(y=`50%`,x=day,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=stage),alpha=0.2) +
geom_line()

ggsave("~/ci.png",g1,dpi=400,scale=1.5,width=25,height=15,unit="cm")