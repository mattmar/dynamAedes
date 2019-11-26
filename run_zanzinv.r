## Load time series of daily average temperature and distance matrices
setwd("~/GitHub/euaeae")
w <- read.csv('data/t_data250m_res100m_clipped.csv')
ww <- w[,-c(1)] #remove index
asd <- apply(as.matrix(dist(ww[,c("x","y")], "maximum", diag=TRUE, upper=TRUE)),2,function(x) round(x/100,1)*100) #Trasform to meters
ww <- ww[,-c(1:2)] #remove lat long
#pld <- readRDS('data/ld_matrix.RDS')
pld <- read.csv('data/t_dist_250m_res100m.txt')
row.names(pld) <- which(w[,1]%in%pld[,1])
colnames(pld) <- row.names(pld)
pld <- apply(pld,2,function(x) round(x/10000,1)*10000) #make round 100's 
gc()

# Run simulations
source("zanzinv.r")
bla <- zanzinv(temp.matrix=ww, cell.dist.matrix=asd, road.dist.matrix=pld,startd=150,endd=360,n.clusters=2,cluster.type="SOCK",iter=2,intro.adults=100)

### Plot results ###
#trend in time, just one iteration
outlp<-mclapply(bla[[1]], function(x) apply(x, MARGIN=c(1, 2), sum),mc.cores=4)
outlt<-melt(t(sapply(outlp,rowSums)))
outlt$day<-as.integer(row.names(outlt))
ggplot(outlt, aes(x=Var1,y=value,col=as.factor(Var2))) + geom_line()

#trend in space, just one iteration
library(rgdal)
wadults<-cbind(w[2:3],sapply(outlp, function(x) x[3,]))
wlarvae<-cbind(w[2:3],sapply(outlp, function(x) x[2,]))
weggs<-cbind(w[2:3],sapply(outlp, function(x) x[1,]))
roads<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/genova/genova_road3035.shp')
area<-readOGR('/home/matteo/own_data/PoD/topics/aegypti_eu/spatial_data/genova/genova3035.shp')

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