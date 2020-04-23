#Sensitivity analysis
#Derive proportion of successful introduction per scenario
setwd('~/euaeae/sa/prop/')
pofe <- mclapply(list.files(pattern="sa_alg*"), function(x) {

	y <- readRDS(x)
	days <- max(sapply(y,length))
	pe <- length(which(sapply(y,function(x) sum(x[days][[1]]))>0)) / length(sapply(y,function(x) sum(x[days][[1]])))
	gc()
	return(pe)

},mc.cores=3)

si_es<-cbind.data.frame(si=unlist(pofe),sc=gsub("[A-Z|\\.S]","",gsub(".*_phe","",list.files(pattern="sa_alg*"))))
si_es<-si_es[order(si_es$sc),]

#Derive maximum median abundance per scenario
setwd('~/euaeae/sa/prop/')
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

###Derive the maximum density and invaded area per scenario
library(parallel)
library(ggplot2)
str=105
endr=1037#732#
days=932

# Abundance
cip <- mclapply(list.files(pattern="sa_alg*"),function(x) {
	y<-readRDS(x);
	z<-rbind.data.frame(returnci(y,st=1,days=days),returnci(y,st=2,days=days),returnci(y,st=3,days=days));
	gc()
	return(z)
},mc.cores=2)

# Number of cells
cic <- mclapply(list.files(pattern="sa_alg*"),function(x) {
	y<-readRDS(x);
	z<-rbind.data.frame(returnnc(y,days=days));
	gc()
	return(z)
},mc.cores=2)

#Reshape abundance and add cells
afi<-lapply(cip, function(x) {x$stage<-as.factor(x$stage); levels(x$stage)<-c("Egg","Immature","Adult"); return(x)})
afi<-lapply(afi, function(x) {x$day<-as.Date(x$day,origin="2017-05-15"); return(x)})
afi<-lapply(1:11, function(x) {afi[[x]]$sc<-rep(gsub("[A-Z|\\.S]","",gsub(".*_pda","",list.files(pattern="sa_alg*")[x])),932*3); return(afi[[x]])})
aci<-lapply(cic, function(x) rbind.data.frame(x,x,x))
afi<-do.call(rbind.data.frame,afi)
aci<-do.call(rbind.data.frame,aci)
names(aci) <- c("25%c","50%c","75%c","day")
afi<-cbind.data.frame(afi,aci[,1:3])
afi$year<-format(afi$day,"%Y")
afi$d50<-afi$`50%`/(afi$`50%c`*6.25)
afi$d25<-afi$`25%`/(afi$`25%c`*6.25)
afi$d75<-afi$`75%`/(afi$`75%c`*6.25)

#Check day of peak adult abundance per hectar per scenario
da_sa <- aggregate(.~sc,afi[afi$stage%in%"Adult",c(6,11:13)],"max")
ci_sa <- aggregate(.~sc,afi[afi$stage%in%"Adult",c(7:9,6)],"max")

#########
#Derive distance betwen point of introduction and cells with at least one individual-do this on the cluster.
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
cc <- mclapply("algeciras_temps11-19_corr.RDS", function(x) readRDS(x)[,2:3],mc.cores=5) #coordinates
ccc <- lapply(cc, function(x) as.data.frame(coordinates(spTransform(SpatialPoints(coords = x, proj4string = CRS("+init=epsg:3035")), CRS("+init=epsg:4326")))))
ccc <- lapply(1,function(x) {ccc[[x]]$distm<-as.integer(distm(introcell[x,1:2],ccc[[x]]));return(ccc[[x]])})

#Function to save only cells with at least 1 mosquito and add distance from introduction cell
setwd('~/euaeae/sa/prop/')
toload<-list.files(pattern="*RDS")
zanzred <- mclapply(1:6, function(x) {
	lcip <-readRDS(toload[x])
	lapply(lcip, function(d) {
		sapply(d, function(y) {
			z <- data.frame(y[,which(colSums(y)>0)]); 
			colnames(z) <- which(colSums(y)>0)
			z[4,]<-ccc[[1]]$distm[which(colSums(y)>0)];
			cat(paste("\n object ",x," is done\n"))
			return(z)
			gc()
		})
	})
},mc.cores=1)

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

asp <- mclapply(zanzred, function(x) cbind.data.frame(returndis(x,days=932),day=1:932),mc.cores=3)
asp <- lapply(asp, function(x) {x$day<-as.Date(x$day,origin="2017-05-15"); return(x)})
asp <- lapply(1:6, function(x) {asp[[x]]$sc<-rep(gsub("[A-Z|\\.S]","",gsub(".*_i200_e","",list.files(pattern="sa_alg*")[x])),932); return(asp[[x]])})
asp <- do.call(rbind.data.frame,asp)

asp[,1:3]<-asp[,1:3]-856

dd_sa <- aggregate(.~sc,asp[,c(1:3,5)],"max")

#Merge outputs and save
outdf <- cbind.data.frame(si_es,da_sa,ci_sa,dd_sa)
write.csv(outdf,"~/outdf_sa_pda.csv")

#Plot
#Egg survival
es<-read.csv('/home/matteo/own_data/PoD/topics/aegypti_eu/output/sa/csv/outdf_sa_es.csv')[,c(3,2,5,6,7,9,10,11,13,14,15)]
names(es) <- c("sc","is","d50","d25","d75","c25","c50","c75","dd25","dd50","dd75")
es$sc<-c(0.9,seq(0.91,0.99,0.01),1.0)
es[,c("c50","c25","c75")]<-es[,c("c50","c25","c75")]*6.25

esm<-cbind.data.frame(melt(es,id.vars="sc",measure.vars=c("d25","c25","dd25","is")),
	v25=melt(es,id.vars="sc",measure.vars=c("d25","c25","dd25","is"),value.name="v25")[,c(3)],
	v50=melt(es,id.vars="sc",measure.vars=c("d50","c50","dd50","is"),value.name="v50")[,c(3)],
	v75=melt(es,id.vars="sc",measure.vars=c("d75","c75","dd75","is"),value.name="v75")[,c(3)]
	)
esm[which(esm$variable%in%"is"),c("v25","v75")]<-NA

levels(esm$variable) <- c("Maximum\ndensity (female/ha)", "Maximum\ninvaded area (ha)", "Maximum\ndispersal distance (m)", "Probability\nof successful invasion")

esm$variable<-factor(esm$variable,levels=levels(esm$variable)[c(4,1,2,3)],ordered=T)
esm$g <- "Daily p of egg survival"

ggplot(esm) +
geom_point(aes(x=sc,y=v50)) +
geom_line(aes(x=sc,y=v50)) +
geom_errorbar(aes(x=sc,ymin=v25,ymax=v75),col="grey") +
facet_wrap(~variable,scales="free_y") +
ylab("") +
xlab("Daily probability of egg survival") +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig02_sm.pdf",dpi=600,scale=2,width=11,height=11,unit="cm")


#probability of passive dispersal
pda<-read.csv('/home/matteo/own_data/PoD/topics/aegypti_eu/output/sa/csv/outdf_sa_pda.csv')[,c(3,2,5,6,7,9,10,11,13,14,15)]
pda$sc <- seq(0.0018,0.0108,length=11)[c(1:8,10,9,11)]
names(pda) <- c("sc","is","d50","d25","d75","c25","c50","c75","dd25","dd50","dd75")

pda[,c("c50","c25","c75")]<-pda[,c("c50","c25","c75")]*6.25

pdam<-cbind.data.frame(melt(pda,id.vars="sc",measure.vars=c("d25","c25","dd25","is"),value.name="v25"),
	v25=melt(pda,id.vars="sc",measure.vars=c("d25","c25","dd25","is"),value.name="v25")[,c(3)],
	v50=melt(pda,id.vars="sc",measure.vars=c("d50","c50","dd50","is"),value.name="v50")[,c(3)],
	v75=melt(pda,id.vars="sc",measure.vars=c("d75","c75","dd75","is"),value.name="v75")[,c(3)]
	)
pdam[which(pdam$variable%in%"is"),c("v25","v75")]<-NA

levels(pdam$variable) <- c("Maximum\ndensity (female/ha)", "Maximum\ninvaded area (ha)", "Maximum\ndispersal distance (m)", "Probability\nof successful invasion")

pdam$variable<-factor(pdam$variable,levels=levels(pdam$variable)[c(4,1,2,3)],ordered=T)
pdam$g <- "Daily p of passive dispersal"

ggplot(pdam) +
geom_point(aes(x=sc,y=v50)) +
geom_line(aes(x=sc,y=v50)) +
geom_errorbar(aes(x=sc,ymin=v25,ymax=v75),col="grey") +
facet_wrap(~variable,scales="free_y") +
ylab("") +
xlab("Probability of passive dispersal") +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig02_sm.pdf",dpi=600,scale=2,width=11,height=11,unit="cm")

#probability of hatching
pha<-read.csv('/home/matteo/own_data/PoD/topics/aegypti_eu/output/sa/csv/outdf_sa_phe.csv')[,c(3,2,5,6,7,9,10,11,13,14,15)]
pha$sc <- as.numeric(paste(0,sapply(strsplit(as.character(pha$sc),split="_",fixed=TRUE), function(x) {x[[2]]}),sep="."))*10
#pha[2,9:11]<-c(400.01,998.67,1348)
names(pha) <- c("sc","is","d50","d25","d75","c25","c50","c75","dd25","dd50","dd75")

pha[,c("c50","c25","c75")]<-pha[,c("c50","c25","c75")]*6.25

pham<-cbind.data.frame(melt(pha,id.vars="sc",measure.vars=c("d25","c25","dd25","is"),value.name="v25"),
	v25=melt(pha,id.vars="sc",measure.vars=c("d25","c25","dd25","is"),value.name="v25")[,c(3)],
	v50=melt(pha,id.vars="sc",measure.vars=c("d50","c50","dd50","is"),value.name="v50")[,c(3)],
	v75=melt(pha,id.vars="sc",measure.vars=c("d75","c75","dd75","is"),value.name="v75")[,c(3)]
	)
pham[which(pham$variable%in%"is"),c("v25","v75")]<-NA

levels(pham$variable) <- c("Maximum\ndensity (female/ha)", "Maximum\ninvaded area (ha)", "Maximum\ndispersal distance (m)", "Probability\nof successful invasion")

pham$variable<-factor(pdam$variable,levels=levels(pdam$variable)[c(4,1,2,3)],ordered=T)
pham$g <- "Daily p of egg hatching"

ggplot(pham) +
geom_point(aes(x=sc,y=v50)) +
geom_line(aes(x=sc,y=v50)) +
geom_errorbar(aes(x=sc,ymin=v25,ymax=v75),col="grey") +
facet_wrap(~variable,scales="free_y") +
ylab("") +
xlab("Probability of passive dispersal") +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig02_sm.pdf",dpi=600,scale=2,width=11,height=11,unit="cm")

#Introduced eggs
prop<-read.csv('/home/matteo/own_data/PoD/topics/aegypti_eu/output/sa/csv/outdf_sa_prop.csv')
prop$sc <- as.numeric(prop$sc)
#pha[2,9:11]<-c(400.01,998.67,1348)
names(prop) <- c("si","sc","d50","d25","d75","c25","c50","c75","dd25","dd50","dd75")

prop[,c("c50","c25","c75")]<-prop[,c("c50","c25","c75")]*6.25

propm<-cbind.data.frame(melt(prop,id.vars="sc",measure.vars=c("d25","c25","dd25","si"),value.name="v25"),
	v25=melt(prop,id.vars="sc",measure.vars=c("d25","c25","dd25","si"),value.name="v25")[,c(3)],
	v50=melt(prop,id.vars="sc",measure.vars=c("d50","c50","dd50","si"),value.name="v50")[,c(3)],
	v75=melt(prop,id.vars="sc",measure.vars=c("d75","c75","dd75","si"),value.name="v75")[,c(3)]
	)
propm[which(propm$variable%in%"si"),c("v25","v75")]<-NA

levels(propm$variable) <- c("Maximum\ndensity (female/ha)", "Maximum\ninvaded area (ha)", "Maximum\ndispersal distance (m)", "Probability\nof successful invasion")

propm$variable<-factor(propm$variable,levels=levels(propm$variable)[c(4,1,2,3)],ordered=T)
propm$g <- "Number of introduced eggs"

ggplot(propm) +
geom_point(aes(x=sc,y=v50)) +
geom_line(aes(x=sc,y=v50)) +
geom_errorbar(aes(x=sc,ymin=v25,ymax=v75),col="grey") +
facet_wrap(~variable,scales="free_y") +
ylab("") +
xlab("Probability of passive dispersal") +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig02_sm.pdf",dpi=600,scale=2,width=11,height=11,unit="cm")


#Plot all together
outam<-rbind.data.frame(
	cbind.data.frame(esm,vl=rep(0.99,nrow(esm))),
	cbind.data.frame(pdam,vl=rep(0.0051,nrow(pdam))),
	cbind.data.frame(pham,vl=rep(0.076,nrow(pham))),
	cbind.data.frame(propm,vl=rep(200,nrow(prop)))
	)


outam$g<-as.factor(outam$g)
outam$g<-factor(outam$g,levels=levels(outam$g)[c(4,2,1,3)],ordered=T)


ggplot(outam) +
geom_point(aes(x=sc,y=v50)) +
geom_line(aes(x=sc,y=v50)) +
geom_vline(aes(xintercept=vl),col="red",linetype=2) +
geom_errorbar(aes(x=sc,ymin=v25,ymax=v75),col="grey") +
facet_grid(variable~g,scales="free") +
labs(x=NULL, y=NULL) +
theme_bw() +
theme(axis.text.x = element_text(angle=0,size=12), axis.text.y = element_text(size=12),axis.title.y = element_text(size=12), legend.position="top",legend.direction="horizontal",legend.text = element_text(size=12),legend.title = element_text(size=1,face="bold"),strip.text.x = element_text(size = 16),strip.text.y = element_text(size = 16))

ggsave("~/own_data/PoD/topics/aegypti_eu/manuout/fig02_sm.pdf",dpi=600,scale=1.4,width=29,height=21,unit="cm")
