## Demographic model for an organism with three life stages (egg, immature [instars and pupa] and adult).
# Spatially explicit; T-dependent daily survival and development rates.
#.a=adult, i.=immature, e.=egg
#.p=probability, .r=rate, .n=number, .v=vector, .m=matrix, .a=array, .f=function

### Temporary lines (will be removed later) ###
## Simulate/load time series of daily average temperature
setwd("~/GitHub/euaeae")
library(actuar)
#library(VFS)
#data("weather") 
#weather.param.p <- wth.param(weather, method = "poisson", llim = 0.3)
#temp1 <- temperature(365*1, weather.param.p)[90:365]
#temp1<-rep(15,365)
#temp1<-read.csv("t_timeseries.csv")$temp[140:363]
# Grid for spatial movement (200*50 correspond to a 10kmx10km extent)
#w <- expand.grid(latcoords = seq(from = 0, by = 200, l = 50),lngcoords = seq(from = 0, by = 200, l = 50))
w <- read.csv('/home/matteo/GitHub/euaeae/t_data.csv')
ww <- w[,-c(1)]
names(ww)[c(1,2)] <- c("lon","lat")
k<-expand.grid(0,c(0.015,-0.015),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
www<-ddply(ww,.(lon,lat),function(x)as.list(x)+k)
www<-www[-c(1:2)]
plot(www[,1]*52000,www[,2]*52000)
points(ww[,1]*52000,ww[,2]*52000,col="red")
names(www)[c(1,2)] <- c("lon","lat")

#Make a matrix of distances
d <- round(as.matrix(dist(www[,c("lon","lat")], "maximum", diag=TRUE, upper=TRUE)*52000/100000),3)*100000 #Trasform from degress to meters?

### Preamble: define variable for the model ###
## Time
days <- ncol(www)-2 #n of day to simulate
## Space
p.life.a <- array(0,c(3,nrow(www),5))
space <- nrow(www)
## Vector of propagules to initiate the life cycle
e.intro.n <- rep(0,space); e.intro.n[sample(1:space,1)] <-100 #n of introduced eggs
i.intro.n <-  0 #n of introduced immatures
a.intro.n <-  0 #n of introduced adults
# Matrix for immature sub-compartments
i.surv.m  <- matrix(0,ncol=5,nrow=2) 
outl <- list(NA)
i.temp.v <- 0

### Life cycle ###
for (day in 1:days) {

### Header: load temperature dependent functions ###
## Gonotrophic cycle
	source("./lc/a.gono_rate.f.r")
## Derive daily rate for gonotrophic cycle,i.e. blood meal to oviposition. 
	a.gono.r <- a.gono_rate.f(www[,day])
## Transform rate in daily probabiltiy to terminate the gonotrophic cycle.
	a.gono.p <- 1-exp(-(a.gono.r))
## Oviposition rate
	source("./lc/a.ovi_rate.f.r")
## Derive oviposition rate, i.e., number of eggs laid per female per day 
	a.batc.n <- a.ovi_rate.f(www[,day])
## Adult survival
	source("./lc/a.surv_rate.f.r")
## Derive daily adult female survival rate
	a.surv.r <- a.surv_rate.f(www[,day])
## Transform rate in daily probabiltiy to survive.
	a.surv.p <- exp(-a.surv.r)
## Add difference between lab and field survival (from Brady et al. 2014)
	a.surv.p <- ifelse(a.surv.p<0,0,a.surv.p)
## Immature survival
	source("./lc/i.surv_rate.f.r")
## Derive daily immature survival rate
	i.surv.r <- i.surv_rate.f(www[,day])
## Transform rate in daily probabiltiy to survive.
	i.surv.p <- exp(-i.surv.r)
## Immature emergence
	source("./lc/i.emer_rate.f.r")
## Derive daily immature emergence rate
	i.emer.r <- i.emer_rate.f(www[,day])
## Transform rate in daily probabiltiy to survive.
	i.emer.p <- exp(-i.emer.r)
## Probability of dispersal
	f.disp.p<-dlnorm(seq(0,600,100),meanlog=4.95,sdlog=0.66)
## Probability of egg survival
	e.surv.p <- 0.99
## Probability of hatching
	e.hatc.p <- exp(-rgamma(length(www[,day]),shape=(7.6/7.4)^2,rate=(7.6/7.4^2))) #egg hatching probability (Soares-Pinheiro et al. 2015)

## Egg compartment ##
# Fixed probabilities
	p.life.a[1,,1] <- if(day==1) e.intro.n else p.life.a[1,,1]+a.egg.n
	
# Random binomial draw to find numbers of eggs that die, survive or hatch
	e.fate.v <- sapply(1:space, function(x){rmultinom(1,p.life.a[1,x,1],c((1-e.surv.p),e.surv.p, e.hatc.p[x]))})
	p.life.a[1,,1] <- e.fate.v[2,] #surviving eggs which did not hatch
	e.hatc.n <- e.fate.v[3,] #eggs hatched
	
## Immature compartment ##
# i.surv.m has five sub-compartments representing days from hatching; an immature can survive/die for the first four days after hatching, from the fifth day on, it can survive/die/emerge. 
	p.life.a[2,,5] <- if(day==1) i.intro.n else 0
	p.life.a[2,,2:5] <- t(sapply(1:space, function(x) {sapply(if(day==1) 0 else p.life.a[2,x,1:4],rbinom,n=1,p=i.surv.p[x])})) #updated n of immatures older than 1day
	p.life.a[2,,1] <- e.hatc.n #Add new immatures 1d old
	p.life.a[2,,5] <- p.life.a[2,,5] + i.temp.v
# Random binomial draw to find numbers of immature 5d+ old that not emerge or emerge
	i.fate.v <- sapply(1:space, function(x){rmultinom(1,p.life.a[2,x,5], c(i.surv.p[day], i.emer.p[x]))})
	i.temp.v <- i.fate.v[1,] #immature surviving not emerging
	i.emer.n <- i.fate.v[2,] #immature emerging

## Adult compartment ##
# a.surv.v has five sub-compartments representing: adults in day 1,2,3 of oviposition; 2d+ old adults host-seeking and non ovipositing; 1d old adults, non-laying and non-dispersing. 
	p.life.a[3,,5] <- if(day==1) a.intro.n else 0
# Remove males adult from the newly emerged adults
	p.life.a[3,,5] <- sapply(1:space,function(x) rbinom(1,i.emer.n[x], 0.5))
# Find number of females which pass from host-seeking to ovipositing
	p.life.a[3,,1] <- sapply(1:space, function(x) rbinom(1, p.life.a[3,x,4], a.gono.p[x])) + p.life.a[3,,1]
	p.life.a[3,,4] <- p.life.a[3,,4] - p.life.a[3,,1]
# Find number of eggs laid by ovipositing females
	a.egg.n <- sapply(1:space, function(x) sum(rpois(p.life.a[3,x,1:3][which(p.life.a[3,x,1:3]>0)], a.batc.n[x])))
# Find number of adult females surviving
	p.life.a[3,,1:4] <- t(sapply(1:space, function(x) {sapply(if(day==1) 0 else p.life.a[3,x,1:4],rbinom,n=1,p=a.surv.p[x])}))
# Find how many host-seeking females move and where
	f.disp.v <- sapply(1:space, function(x) rmultinom(1,p.life.a[3,x,4],cbind(f.disp.p)))
# Find cells which have at least 1 host-seeking dispersing adult (cells of origin)
	f.ocel.v <- which(colSums(f.disp.v)>0)
# Find a landing cell for each distance at which a set of adult females disperse 
	for (i in f.ocel.v) {
		a.plan.l <- sapply(which(which(f.disp.v[,i]>0)*100<max(d[i,]))*100, function(x) {dx<-d[i,which(d[i,]==ifelse(x==200,100,x))]})#Find set of dispersing cells for each cell of origin
		if(length(unlist(a.plan.l))>0) a.land.v <- as.integer(names(sapply(a.plan.l, sample,1,replace=TRUE))) else next() #Randomly choose a cell of landing for each dispersing distance
		p.life.a[3,a.land.v,4] <- f.disp.v[which(which(f.disp.v[,i]>0)*100<max(d[i,])),i] #Add dispersing individuals in the landing cell
	}
# Make a new host-seeking compartment and slide ovipositing female status for new day
	p.life.a[3,,4] <- p.life.a[3,,3] + p.life.a[3,,4] + p.life.a[3,,5]
	p.life.a[3,,2:3] <- p.life.a[3,,1:2]
	p.life.a[3,,1] <- 0
	outl[[day]] <- p.life.a
	message("day",day,"\n")
# If the sum of all compartments and cells is 0, the population is exinct, stop the life cycle.
	if(sum(p.life.a)==0) stop("Exinct!")
}

### Plot results ###
#trend in time
outlp<-lapply(outl, function(x) apply(x, MARGIN=c(1, 2), sum))
outlt<-melt(t(sapply(outlp,rowSums)))
outlt$day<-as.integer(row.names(outlt))
ggplot(outlt, aes(x=Var1,y=value,col=as.factor(Var2))) + geom_line()

#trend in space
wadults<-cbind(www[1:2]*6556.751,sapply(outlp, function(x) x[3,]))

for (i in 1:(ncol(wadults)-2)) {
	plot(wadults[,1],wadults[,2],cex=0.1,col="green")
	text(x=900,y=1000,as.Date(i+130, origin = "2017-01-01"),col="red")
	points(wadults[,1],wadults[,2],cex=wadults[,i+2]/5,fill="red")
	Sys.sleep(0.1)
}