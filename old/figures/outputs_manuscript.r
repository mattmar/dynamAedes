#Question #1: could Aedes aegypti introduced via ship trade in the major European harbours establish viable populations
##Find most favourable period of the year for a successfull introduction
setwd('/home/matteo/own_data/PoD/topics/aegypti_eu/output/preiter')
iterciall<-do.call(rbind.data.frame,lapply(list.files(pattern="*g1.RDS"),readRDS))
iterciall$harbour<-as.factor(rep(c("Algeciras","Barcelona","Genova","Rotterdam","Venezia"), each=nrow(iterciall)/5))
iterciall$harbour<-factor(iterciall$harbour,levels=levels(iterciall$harbour)[c(1,2,3,5,4)],ordered=TRUE)
iterciall$month<-factor(month.abb[iterciall$month],levels=month.abb)

#Get the number of adults at peak
iterciall[iterciall$stage%in%"Adult",][which.max(iterciall[iterciall$stage%in%"Adult",]$`50%`),]

ggplot(iterciall, aes(y=`50%`+1,x=date,group=stage,col=stage)) + 
#geom_ribbon(aes(ymin=`2.5%`+0.01,ymax=`97.5%`+0.01,fill=stage),alpha=0.1) +
geom_line() +
facet_grid(month~harbour) +
labs(x="Month from introduction",y="Abundance (Log10)") +
scale_x_date(date_breaks = "2 months", date_labels = "%b",limits=c(as.Date("2017-01-01"),as.Date("2018-01-01"))) +
theme(axis.text.x = element_text(angle = 45,size=8),axis.text.y = element_text(size=8), legend.position="bottom",legend.text = element_text(size=15),legend.title = element_text(size=12)) +
scale_y_continuous(labels=function(x) {format(x, scientific = TRUE)},trans="log10")
ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/intromonth_new.png",dpi=400,scale=1,width=21,height=29,unit="cm")

#Get max per harbour
aggregate(`97.5%`~month,data=iterciall[iterciall$stage%in%"Adult"&iterciall$harbour%in%"Barcelona",],"max")

#Get the number of adults during the last day 2017-12-21
ld<-iterciall[iterciall$date%in%as.Date("2017-12-21"),]
aggregate(`50%`~month+harbour,ld,"sum",na.pass=T)

#This is the probability that an introduction is successfull at the end of the simulated period 
##Find the probability od successful introduction
setwd('~/euaeae/sa/prop/')

# Apply to RDS
pofe <- mclapply(list.files(pattern="sa_bar*"), function(x) {

	y <- readRDS(x)
	days <- max(sapply(y,length))
	pe <- length(which(sapply(y,function(x) sum(x[days][[1]]))>0)) / length(sapply(y,function(x) sum(x[days][[1]])))
	gc()
	return(pe)

},mc.cores=1)

cbind.data.frame(unlist(pofe),list.files(pattern="sa_bar*"))

#Question #3
##Need to introduce different numbers of propagules

#Question #4: What is the space-time trend which Ae. aegypti populations will follow once established
##Time trend
###Derive the maximum number of day in a simulation
library(parallel)
library(ggplot2)
str=105
endr=1037#732#
days=932

setwd('~/euaeae/fulliter/')

#Function to derive CI of abundance per cell per day
returnci<-function(outl=NA,st=1,cores=5,days=0,ppp=c(0.25,0.50,0.75)){
	out<-apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=ppp,na.rm=T);
	colnames(out)<-NULL
	outo<-rbind.data.frame(out,
		stage=rep(st,nrow(out)),
		day=as.factor(1:ncol(out)))
	return(t(outo))
}

#Function to derive CI of number of cell invaded per day
returnnc<-function(outl=NA,cores=5,days=0,ppp=c(0.25,0.50,0.75)){
	out<-apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
		lapply(1:days, function(y) {
			if(y<=length(x)) {length(which(x[[y]]>0))} else {NA}})},mc.cores=cores)),2,quantile,probs=ppp,na.rm=T);
	colnames(out)<-NULL
	outo<-rbind.data.frame(out,
		day=as.factor(1:ncol(out)))
	return(t(outo))
}

#Derive CI of abundance and invaded cells per stage - do this on the cluster
# Abundance
cip <- mclapply(list.files(pattern="*RDS"),function(x) {
	y<-readRDS(x);
	z<-rbind.data.frame(returnci(y,st=1,days=days),returnci(y,st=2,days=days),returnci(y,st=3,days=days));
	gc()
	return(z)
},mc.cores=3)

# Number of cells
cic <- mclapply(list.files(pattern="*RDS"),function(x) {
	y<-readRDS(x);
	z<-rbind.data.frame(returnnc(y,days=days));
	gc()
	return(z)
},mc.cores=3)

saveRDS(cip,"~/allfulliter100e.RDS")
saveRDS(cic,"~/allfulliter100e_cells.RDS")

# laptop: load time trend and reshape
afi<-readRDS('/home/matteo/own_data/PoD/topics/aegypti_eu/output/fulliter/allfulliter100e.RDS')
aci<-readRDS('/home/matteo/own_data/PoD/topics/aegypti_eu/output/fulliter/allfulliter100e_cells_new.RDS')

afi<-lapply(afi, function(x) {x$stage<-as.factor(x$stage); levels(x$stage)<-c("Egg","Immature","Adult"); return(x)})
#aci<-lapply(aci, function(x) {x$stage<-as.factor(x$stage); levels(x$stage)<-c("Egg","Immature","Adult"); return(x)})

afi[c(1,3,5)]<-lapply(afi[c(1,3,5)], function(x) {x$day<-as.Date(x$day,origin="2017-05-15"); return(x)})
afi[2]<-lapply(afi[c(2)], function(x) {x$day<-as.Date(x$day,origin="2017-04-15"); return(x)})
afi[4]<-lapply(afi[c(4)], function(x) {x$day<-as.Date(x$day,origin="2017-03-15"); return(x)})

aci<-lapply(aci,function(x) rbind.data.frame(x,x,x))
afi<-do.call(rbind.data.frame,afi)
aci<-do.call(rbind.data.frame,aci)

names(aci) <- c("25%c","50%c","75%c","day")

afi<-cbind.data.frame(afi,aci[,1:3])


#Plot all trends together
plot1 <- 
ggplot(afi, aes(y=(`50%`/(`50%c`*6.25))+1,x=day,group=stage,col=stage)) + 
facet_wrap(~port,ncol=2) +
geom_ribbon(aes(ymin=(`25%`/(`50%c`*6.25))+1,ymax=(`75%`/(`50%c`*6.25))+1,fill=stage),colour="transparent",alpha=0.4,show.legend=FALSE) +
geom_line(linetype=1,size=1) +
scale_color_brewer(palette = "Dark2",direction=-1) +
scale_fill_brewer(palette = "Dark2",direction=-1) +
facet_wrap(~port,ncol=1,scale="free_y") +
labs(x=NULL,y=expression(paste("Abundance ", (ha^-1)))) +
guides(colour=guide_legend(title="")) +
scale_x_date(date_breaks = "3 months", date_labels = "%b\n%y",limits=c(as.Date("2017-03-15"),as.Date("2019-09-30")),expand=c(0.01,0.01)) +
scale_y_continuous(trans="log10") +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16)); plot1

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig01_sm.pdf",plot1,dpi=600,scale=1,width=21,height=27,unit="cm")


#Import temperature to add to the time trend figure
algt <- as.numeric(readRDS('~/GitHub/euaeae/data/algeciras_temps11-19_corr.RDS')[,2193:3229][15664,166:1037])/1000
bart <- as.numeric(readRDS('~/GitHub/euaeae/data/barcelona_temps17-19_corr.RDS')[70616,136:1037])/1000
gent <- as.numeric(readRDS('~/GitHub/euaeae/data/genova_temps17-19_corr.RDS')[46562,166:1037])/1000
rott <- as.numeric(readRDS('~/GitHub/euaeae/data/rotterdam_temps17-19_corr.RDS')[37776,106:1037])/1000
vent <- as.numeric(readRDS('~/GitHub/euaeae/data/venezia_temps17-19_corr.RDS')[27475,166:1037])/1000

#Add temperature of a pixel in the harbour to the main CI dataframe
afi$temp<-c(rep(c(algt,rep(NA,60)),3),rep(c(bart,rep(NA,30)),3),rep(c(gent,rep(NA,60)),3),rep(rott,3),rep(c(algt,rep(NA,60)),3))
afi$port<-as.factor(rep(c("Algeciras","Barcelona","Genova","Rotterdam","Venezia"),each=932*3))
afi$port<-factor(afi$port,levels=levels(afi$port)[c(1,2,3,5,4)],ordered=T)
afi$year<-format(afi$day,"%Y")

#Check day of peak adult abundance per hectar per port
afi$density<-afi$`50%`/(afi$`50%c`*6.25)

maxv<-aggregate(`75%`~stage+port+year,afi[afi$stage%in%"Egg"&afi$port%in%"Barcelona",],"min")

afi[which(afi$`50%`%in%maxv&afi$port%in%"Algeciras"),]


#Space trend
#Derive ditance betwen each point of introduction and each cell which has at least one individual- do this on the cluster
library(geosphere)
library(parallel)
library(sp)

introcell <- coordinates(spTransform(SpatialPoints(matrix(c(
	2928516.524,1592445.244,#alg
	3662783.371,2061013.710,#bar
	4234031.417,2366471.634,#gen
	3911849.229,3223083.767,#rot
	4497026.260,2484877.842 #ven
	),ncol=2,byrow=T), proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))

#Calculate distance between coords of introduction and all other cells
setwd('~/GitHub/euaeae/data/')
cc <- mclapply(list.files(pattern="*RDS"), function(x) readRDS(x)[,2:3],mc.cores=5) #coordinates
ccc <- lapply(cc, function(x) as.data.frame(coordinates(spTransform(SpatialPoints(coords = x, proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))))
ccc <- lapply(1:5,function(x) {ccc[[x]]$distm<-as.integer(distm(introcell[x,1:2],ccc[[x]]));return(ccc[[x]])})

#Function to save only cells with at least 1 mosquito and add distance from introduction cell
setwd('~/euaeae/sa/p')
toload<-list.files(pattern="*RDS")
zanzred <- mclapply(1:length(toload), function(x) {
	lcip <-readRDS(toload[x])
	lapply(lcip, function(d) {
		sapply(d, function(y) {
			z <- data.frame(y[,which(colSums(y)>0)]); 
			colnames(z) <- which(colSums(y)>0)
			z[4,]<-ccc[[x]]$distm[which(colSums(y)>0)];
			cat(paste("\n object ",x," is done\n"))
			return(z)
			gc()
		})
	})
},mc.cores=5)

#Function to derive 95% CI across simulations of distance from initial cell of introduction for each day
returndis <- function(outl=NA,days=0){
	out <- do.call(rbind.data.frame, lapply(1:days, 
		function(x) {
			quantile(unlist(sapply(1:length(outl),
				function(y) {mean(if(x<length(outl[[y]])) as.integer(outl[[y]][[x]][4,]) else NA,na.rm=T)})), probs=c(0.25,0.50,0.75), na.rm=T)
		}))
	names(out) <- c("25%","50%","75%")
	return(out)
}

asp <- mclapply(zanzred, function(x) cbind.data.frame(returndis(x,days=932),day=1:932),mc.cores=5)
saveRDS(asp,"~/allfulliter100e_space.RDS")

#Laptop
##Import and reshape
asp <- readRDS('/home/matteo/own_data/PoD/topics/aegypti_eu/output/fulliter/paper/allfulliter100e_space.RDS')
asp[c(1,3,5)]<-lapply(asp[c(1,3,5)], function(x) {x$day<-as.Date(x$day,origin="2017-05-15"); return(x)})
asp[2]<-lapply(asp[c(2)], function(x) {x$day<-as.Date(x$day,origin="2017-04-15"); return(x)})
asp[4]<-lapply(asp[c(4)], function(x) {x$day<-as.Date(x$day,origin="2017-03-15"); return(x)})
asp<-do.call(rbind.data.frame,asp)
asp$port<-as.factor(rep(c("Algeciras","Barcelona","Genova","Rotterdam","Venezia"),each=932))
asp$port<-factor(asp$port,levels=levels(asp$port)[c(1,2,3,5,4)],ordered=T)
asp[asp$port%in%"Algeciras",1:3]<-(asp[asp$port%in%"Algeciras",1:3]-1260)
asp[asp$port%in%"Barcelona",1:3]<-(asp[asp$port%in%"Barcelona",1:3]-2520)
asp[asp$port%in%"Genova",1:3]<-(asp[asp$port%in%"Genova",1:3]-1150)
asp[asp$port%in%"Venezia",1:3]<-(asp[asp$port%in%"Venezia",1:3]-600)
asp[asp$port%in%"Rotterdam",1:3]<-(asp[asp$port%in%"Rotterdam",1:3]-1450)

#Plot trend of population dispersal in space
plot2 <- ggplot(asp, aes(y=`50%`,x=day)) + 
geom_line() +
geom_ribbon(aes(ymin=`25%`,ymax=`75%`),colour="grey90",alpha=0.4) +
scale_color_brewer(palette = "Dark2") +
scale_fill_brewer(palette = "Dark2") +
facet_wrap(~port,ncol=1) +
labs(x=NULL,y="Distance (km) from port of introduction") +
scale_y_continuous(limits=c(0,3800),labels=function(x) format(x/1000, scientific = FALSE)) +
scale_x_date(date_breaks = "3 months", date_labels = "%b\n%y",limits=c(as.Date("2017-03-15"),as.Date("2019-09-30")),expand=c(0.01,0.01)) +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16)); plot2

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig06.pdf",plot2,dpi=600,scale=1,width=15,height=29,unit="cm")

#Plot all trends together
plot3 <- 
ggplot(afi[afi$stage%in%"Adult",], aes(y=(`50%`/(`50%c`*6.25)),x=day,group=stage,col="Adult females")) + 
geom_line(aes(y=temp*4,col="Temperature"),linetype=3) + 
geom_line(aes(y=(`50%c`*6.25)/6,col="Invaded area (ha/5)"),linetype=1) +
facet_wrap(~port,ncol=2) +
geom_ribbon(aes(ymin=(`25%`/(`50%c`*6.25)),ymax=(`75%`/(`50%c`*6.25)),fill=stage),colour="transparent",alpha=0.4,show.legend=FALSE) +
geom_ribbon(aes(ymin=(`25%c`*6.25)/6,ymax=(`75%c`*6.25)/6,fill="black"),colour="transparent",alpha=0.4,show.legend=FALSE) +
geom_line(linetype=1,size=1) +
scale_color_manual(values = c("#1B9E77","black","#A6CEE3","grey60")) +
scale_fill_manual(values = c("#1B9E77","black","#A6CEE3","grey60")) +
facet_wrap(~port,ncol=1,scale="free_y") +
labs(x=NULL,y=expression(paste("Adult female abundance ", (ha^-1)))) +
guides(colour=guide_legend(title="")) +
scale_y_continuous(sec.axis = sec_axis(name="Temperature (°C)",trans= ~ ./4),expand=c(0,0)) +
scale_x_date(date_breaks = "3 months", date_labels = "%b\n%y",limits=c(as.Date("2017-03-15"),as.Date("2019-09-30")),expand=c(0.01,0.01)) +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16)); plot3

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig05.pdf",plot3,dpi=600,scale=1,width=21,height=27,unit="cm")

# New legend position
library(tidyverse)
library(lemon)

shift_legend3 <- function(p) {
	pnls <- cowplot::plot_to_gtable(p) %>% gtable::gtable_filter("panel") %>%
	with(setNames(grobs, layout$name)) %>% purrr::keep(~identical(.x,zeroGrob()))

	if( length(pnls) == 0 ) stop( "No empty facets in the plot" )

		lemon::reposition_legend( p, "center", panel=names(pnls) )
}

# Move the ledeng in the empty facet
plot4<-shift_legend3(plot3)
ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/abunspacetrend_distance.pdf",plot4,dpi=600,scale=1.7,width=21,height=21,unit="cm")

##Map population abundance
library(tidyverse)
setwd("~/GitHub/euaeae/data/")
#Load coordinates for each cell name
portc<-mclapply(list.files(pattern="*corr.RDS"), function(x) {y<-readRDS(x)[,1:3]; return(y)})
portc1<-lapply(portc, function(x) {names(x)[1]<-"cell_id"; return(x)})
names(portc1)<-c("Algeciras","Barcelona","Genova","Rotterdam","Venezia")

#Derive the average number of mosquitoes in each cell across time
#abxcel<-mclapply(zanzred[2:3], function(y){lapply(y, function(x) lapply(x, function(z) {z[4,]<-NA;colSums(z,na.rm=T)} ))},mc.cores=5)
abxcel<-mclapply(zanzred, function(y){lapply(y, function(x) {lapply(x, function(z) {out<-z[3,];return(out)})})},mc.cores=5)
abxcelred<-mclapply(abxcel, function(y) {bind_rows(lapply(y, function(x) aggregate(unlist(x),list(names(unlist(x))),"mean")))},mc.cores=5)
#Derive the interquartile of the average distribution across iterations
statsxcel <- mclapply(abxcelred, function(y) {aggregate(x~Group.1,data=y,"quantile",prob=c(0.25,0.5,0.75))},mc.cores=5)
statsxcel1<-lapply(statsxcel, function(x) {x<-cbind.data.frame(row_id=as.numeric(as.character(x[,1])),unname(x[,2]));return(x)})
statsxcel1<-lapply(1:5, function(x) {statsxcel1[[x]]$cell_id<-portc1[[x]]$cell_id[statsxcel1[[x]]$row_id];return(statsxcel1[[x]])})
names(statsxcel1)<-c("Algeciras","Barcelona","Genova","Rotterdam","Venezia")

portc1 <- bind_rows(portc1,.id="port")
statsxcel1 <- bind_rows(statsxcel1,.id=c("port"))

#Merge coordinates with statistics per cell and beautify it
statsout<-merge(portc1,statsxcel1,by=c("port","cell_id"),all.x=T)
statsout$row_id<-NULL
statsout$port.y<-NULL
names(statsout)[5:7]<-c("CI25%","CI50%","CI75%")
statsout[5:7]<-round(statsout[5:7],1)*2
#statsout[,6]<-ifelse(statsout[,6]<0.5,-999,statsout[,6])
#statsout[,6]<-ifelse(is.na(statsout[,6]),-999,statsout[,6])
#statsout<-cbind.data.frame(statsout,log=round(log10(statsout[5:7]),2))

#Export the file for QGIS
write.csv(statsout,"~/statsout.csv",row.names=FALSE,na="NA")

#For later checks
statsout<-read.csv('/home/matteo/own_data/PoD/topics/aegypti_eu/output/statsout.csv')
#Get max abundance in a cell per port
aggregate(`CI50.`~port,statsout,"max")
#Get number of cells with at least one adult
aggregate(`CI50.`~port,statsout,function(x) length(which(x>0)))

