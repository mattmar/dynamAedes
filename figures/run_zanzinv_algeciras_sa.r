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
introvector <- c(15664,15665,15666,15667,15685,15686,15687,15688,15704,15705,15706,15707,15708,15723,15724,15725,15726,15727,15741,15742,15743,15744,15745,15757,15758,15759,15760,15761,15762,15774,15775,15776,15777,15778,15779,15780,15792,15793,15794,15795,15796,15797,15798,15810,15811,15812,15813,15814,15815,15816,15832,15833,15834,15835,15836,15837,15838,15854,15855,15856,15857,15858,15859,15860,15875,15876,15877,15878,15879,15880,15881,15894,15895,15896,15897,15898,15899,15900,15912,15913,15914,15915,15916,15917,15929,15930,15933,15934)

#Run sensitivity analysis on number of introducted eggs
source("zanzinv_sa.r")
for (egg in 1000) {
	cat("N of introduced eggs: ",egg)
	str=135; endr=1037; ie=egg; nc=25; it=200
	zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld, startd=str,endd=endr,n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=introvector,intro.eggs=ie,compressed.output=TRUE,country="es",e.surv.p=0.99,p_dis_a=0.0051,e_hatch_pa=0.076,suffix=paste("sa_alg_sd",str,"_ed",endr,"_i",it,"_e",ie,sep=""))
}

#Run sensitivity analysis on egg survival
source("zanzinv_sa.r")
for (es in seq(0.90,1.00,length=11)) {
	cat("egg surv p: ",es)
	str=135; endr=1037; ie=100; nc=30; it=50
	zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld, startd=str,endd=endr,n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=introvector,intro.eggs=ie,compressed.output=TRUE,country="es",e.surv.p=es,p_dis_a=0.0051,e_hatch_pa=0.076,suffix=paste("sa_alg_sd",str,"_ed",endr,"_i",it,"_e",ie,"_es_",es,sep=""))
}

#Run sensitivity analysis on the probability that an adult disperses passively
source("zanzinv_sa.r")
for (pda in seq(0.0018,0.0108,length=11)) {
	cat("Adult dispersal p: ",pda,"\n")
	str=135; endr=1037; ie=250; nc=25; it=200
	zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld, startd=str,endd=endr,n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=introvector,intro.eggs=ie,compressed.output=TRUE,country="es",e.surv.p=0.99,p_dis_a=pda,e_hatch_pa=0.076,suffix=paste("sa_alg_sd",str,"_ed",endr,"_i",it,"_e",ie,"_pda_",pda,sep=""))
}

#Run sensitivity analysis on the probability that an egg hatches
source("zanzinv_sa.r")
phe <- seq(0.01,1,length=11)
lt <- seq(0.005,0.5,length=11)
for ( i in 1:11 ) {
	cat("Egg probability of hatching: ",phe[i],"\n")
	str=135; endr=1037; ie=250; nc=25; it=200
	zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld, startd=str,endd=endr,n.clusters=nc, cluster.type="SOCK",iter=it,intro.cells=introvector,intro.eggs=ie,compressed.output=TRUE,country="es",e.surv.p=0.99,p_dis_a=0.0051,e_hatch_pa=phe[i],e_hatch_pam=lt[i],suffix=paste("sa_alg_sd",str,"_ed",endr,"_i",it,"_e",ie,"_phe_",phe[i],"_",lt[i],sep=""))
}

#Single parameters
source("zanzinv_sa.r")
str=135; endr=1037; ie=100; nc=1; it=1
str=135; endr=1037; ie=100; nc=23; it=50
zanzout <- zanzinv(temps.matrix=ww, cells.coords=cc, road.dist.matrix=pld, startd=str,endd=endr, n.clusters=nc, cluster.type="SOCK",iter=it, intro.cells=introvector, intro.eggs=ie, compressed.output=TRUE, country="es", e.surv.p=0.99, p_dis_a=0.0051, e_hatch_pa=0.076, e_hatch_pam=0.023, suffix=paste("out_alg_sd",str,"_ed",endr,"_i",it,"_e",ie,sep=""))

#Derive number of time an iteration is still active at the end of next spring
library(parallel)
library(ggplot2)
setwd('~/euaeae/sa/')

esrds <- mclapply(list.files(pattern="*_es_*"),readRDS,mc.cores=11)
pdards <- mclapply(list.files(pattern="*_pda_*"),readRDS,mc.cores=11)
pherds <- mclapply(list.files(pattern="*_phe_*"),readRDS,mc.cores=11)

days <- max(sapply(pherds, function(x) sapply(x,length)))

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

#Calculate 95% CI adult abundance for each scenario in each day
citoplot<-mclapply(pherds, function(x) {rbind.data.frame(returnci(x,1,days=days),returnci(x,2,days=days),returnci(x,3,days=days))},mc.cores=11)
citoplot<-lapply(citoplot,function(x) {x$stage<-as.factor(x$stage);return(x)})
citoplot<-lapply(citoplot,function(x) {x$day<-as.Date(x$day,origin="2017-04-15");return(x)})
citoplot<-lapply(1:11,function(x) {citoplot[[x]]$value<-as.character(gsub(".RDS","",gsub("^.*_phe_","",list.files(pattern="*_phe_*"))))[x];return(citoplot[[x]])})
citoplot<-do.call(rbind.data.frame,citoplot)

g1 <- ggplot(citoplot, aes(y=`50%`,x=day,group=stage,col=stage)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`,fill=stage),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Populazion size") +
facet_wrap(~value) +
xlim(range(seq(as.Date("2017-04-15"), by = "day", length.out = endr-str))) + 
ylim(0,50000)

ggsave("~/sa_phe_test.png",g1,dpi=400,scale=1.5,width=25,height=15,unit="cm")

###Spatial spread
library(rgdal)
library(geosphere)
w <- readRDS('~/GitHub/euaeae/data/algeciras_temps11-19_corr.RDS')[,c(1:3,2193:3229)] #select years 1=2011;
cc <- w[,c(2:3)] #coordinates
ccc <- coordinates(spTransform(SpatialPoints(coords = cc, proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))
cc$distm <- as.integer(distm(ccc[15664,],ccc))

#save only cells with at least 1 guy
zanzred <- mclapply(pherds, function(y) {lapply(y, function(x) {sapply(x, function(y) {
	z <- data.frame(y[,which(colSums(y)>0)]); 
	colnames(z) <- which(colSums(y)>0)
	z[4,] <- cc$distm[which(colSums(y)>0)]; 
	return(z)
})})},mc.cores=11)

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

plotdistt<-lapply(zanzred, function(x) cbind.data.frame(returndis(x,days=days),day=1:days))
plotdistt<-lapply(plotdistt, function(x) {x$day<-as.Date(x$day,origin="2017-04-15"); return(x)})
plotdistt<-lapply(1:11,function(x) {plotdistt[[x]]$value<-as.numeric(gsub(".RDS","",gsub("^.*_pda_","",list.files(pattern="*_pda_*"))))[x];return(plotdistt[[x]])})
plotdistt<-do.call(rbind.data.frame,plotdistt)

g2 <- ggplot(plotdistt, aes(y=`50%`,x=day)) + 
geom_ribbon(aes(ymin=`2.5%`,ymax=`97.5%`),alpha=0.2) +
geom_line(linetype=1,size=1.5) + labs(y="Distance (m)") +
facet_wrap(~value) +
xlim(range(seq(as.Date("2017-04-15"), by = "day", length.out = endr-str)))

ggsave("~/sa_dt_esa_test.png",g2,dpi=400,scale=1.5,width=25,height=15,unit="cm")