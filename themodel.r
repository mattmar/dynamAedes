## Demographic model for an organism with three life stages (egg, immature [instars and pupa] and adult).
# This simulates a single pixel with n introduced individuals. T-dependent daily survival and development rates.
#.a=adult, i.=immature,e.=egg
#.p=probability, .r=rate, .n=number, .v=vector, .m=matrix, .f=function

### Temporary lines (will be removed later) ###
## Simulate/load time series of daily average temperature
setwd("~/GitHub/euaeae")
#library(VFS)
#data("weather") 
#weather.param.p <- wth.param(weather, method = "poisson", llim = 0.3)
#temp1 <- temperature(365*1, weather.param.p)[90:365]
#temp1<-rep(28,365)
temp1<-read.csv("t_timeseries.csv")$temp[150:363]

### Header: load temperature dependent functions ###
## Gonotrophic cycle
source("./lc/a.gono_rate.f.r")
## Derive daily rate for gonotrophic cycle,i.e. blood meal to oviposition. 
a.gono.r <- a.gono_rate.f(temp1)
## Transform rate in daily probabiltiy to terminate the gonotrophic cycle.
a.gono.p <- 1-exp(-(a.gono.r))
## Oviposition rate
source("./lc/a.ovi_rate.f.r")
## Derive oviposition rate, i.e., number of eggs laid per female per day 
a.batc.n <- ifelse(a.ovi_rate.f(temp1)<0,0,a.ovi_rate.f(temp1))
## Adult survival
source("./lc/a.surv_rate.f.r")
## Derive daily adult female survival rate
a.surv.r <- ifelse(a.surv_rate.f(temp1)<0,0,a.surv_rate.f(temp1))
## Transform rate in daily probabiltiy to survive.
a.surv.p <- exp(-a.surv.r)
## Immature survival
source("./lc/i.surv_rate.f.r")
## Derive daily immature survival rate
i.surv.r <- ifelse(i.surv_rate.f(temp1)<0,0,i.surv_rate.f(temp1))
## Transform rate in daily probabiltiy to survive.
i.surv.p <- exp(-i.surv.r)
## Immature emergence
source("./lc/i.emer_rate.f.r")
## Derive daily immature emergence rate
i.emer.r <- ifelse(i.surv_rate.f(temp1)<0,0,i.surv_rate.f(temp1))
## Transform rate in daily probabiltiy to survive.
i.emer.p <- exp(-i.emer.r)

### Preamble: define variable for the model ###
## Vector of propagules to initiate the life cycle
e.intro.n <- 500 #n of introduced eggs
i.intro.n <-   0  #n of introduced immatures
a.intro.n <-   0  #n of introduced adults
## Time
days <- length(temp1) #n of day to simulate
## Output dataframe where to store numbers of egg, imm and adu per day
res_dt <- as.data.frame(matrix(NA,nrow=days,ncol=3))
names(res_dt) <- c("egg","immature","adult")
#Matrix for immature sub-compartments
i.surv.m  <- matrix(0,ncol=5,nrow=2) 

### Life cycle ###
for (day in 1:days) {

## Egg compartment ##
# Fixed probabilities
	e.surv.p <- 0.99 #egg survival probability
	e.hatc.p <- 0.20 #egg hatching probability
	e.surv.n <- ifelse(day==1,e.intro.n,e.surv.n)
	res_dt[day,1] <- e.surv.n

# Random binomial draw to find numbers of eggs that die, survive or hatch
	e.fate.v <- rmultinom(1, e.surv.n, c((1-e.surv.p), e.surv.p, e.hatc.p))
	e.stay.n <- e.fate.v[2,] #surviving eggs which did not hatch
	e.hatc.n <- e.fate.v[3,] #eggs hatched

## Immature compartment ##
# This compartment has 5 sub-compartments representing days from hatching; an immature can survive/die for the first four days after hatching, from the fifth day on, it can survive/die/emerge. 
	i.stay.n <- ifelse(day==1,i.intro.n,i.stay.n)
	i.surv.m[1,2:5] <- sapply(i.surv.m[2,1:4],rbinom,n=1,p=i.surv.p[day]) #updated n of immatures older than 1day
	i.surv.m[1,5] <- i.surv.m[1,5] + i.stay.n #Add n of immature 5d+ old which did not die the day before (t-1)
	i.surv.m <- apply(i.surv.m, 2, rev) #Reverse the matrix
	i.surv.m[2,1] <- e.hatc.n #Add new immatures 1d old
	res_dt[day,2] <- sum(i.surv.m[2,]) #Save total number of immatures
# Random binomial draw to find numbers of immature 5d+ old that die, survive or emerge
	i.fate.v <- rmultinom(1, i.surv.m[2,5],c((1-i.surv.p[day]), i.surv.p[day], i.emer.p[day]))
	i.stay.n <- i.fate.v[2,] #immature surviving not emerging
	i.emer.n <- i.fate.v[3,] #immature emerging

## Adult compartment ##
	a.surv.n <- ifelse(day==1,a.intro.n,a.surv.n)
# Remove males adult from the population using a random binomial draw	
	a.surv.n <- as.integer(rbinom(1,a.surv.n + i.emer.n,0.5))
# Find how many females are ready to lay eggs using a binomial draw	
	a.laying.n <- rbinom(1,a.surv.n, a.gono.p[day])
# Random draw from binomial distribution to find how many adult female die or survive
	adu_fate.v <- rmultinom(1, a.surv.n, c((1 - a.surv.p[day]), a.surv.p[day]))
	a.surv.n <- adu_fate.v[2,] 
	res_dt[day,3]  <- a.surv.n
# Find number of eggs laid by females by a Poisson draw
	a.egg.n <- rpois(a.laying.n,a.batc.n[day])
	a.egg.n <- ifelse(length(a.egg.n)==0,0,a.egg.n)

# Updated egg reservoir to start the new "day"
	e.surv.n <- e.stay.n + a.egg.n
}

### Plot results ###
res_lf<-cbind.data.frame(melt(cbind.data.frame(res_dt,temp=temp1)),days=rep(1:days,4))
ggplot(res_lf, aes(x=days,y=value,col=variable)) + geom_line()