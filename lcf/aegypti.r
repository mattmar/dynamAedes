#Temperature-dependent life cycle functions for `Aedes aegypti`

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
a.gono_rate.f <- function(temp.new) {
	gono.v <- c(30,15,7,5,4,4,5,5,6,7,35)
	temp.v <- c(15,17,20,22,26,28,31,32,33,33,40)
	lm <- lm(gono.v ~ poly(temp.v,6))
	a.gono.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
	return( a.gono.pred ) 
}
# plot(-15:55,a.gono_rate.f(-15:55),col="red",type="l")
# points(1/gono.v~temp.v,col="blue")

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
a.ovi_rate.f <- function(temp.new) {
	ovi.v <- c(0,0,0.35,1.12,3.37,3.6,6.98,7.6,9.58,7.28,11.22,7.27) * 10
	temp.v <- c(6,10.54,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	lm <- lm(ovi.v ~ poly(temp.v,4))
	a.ovi.pred <- predict(lm,data.frame(temp.v=temp.new))
	a.ovi.pred <- ifelse(a.ovi.pred<0,0,a.ovi.pred)
	return( a.ovi.pred )
}

## Adult female daily mortality rate at different temperatures ##
a.surv_rate.f <- function(temp.new) {
	surv.time <- c(0,0,13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07)
	temp.v <- c(-20,5,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	lm <- lm(surv.time ~ poly(temp.v,3))
	a.surv.rate <- 1-1/predict(lm,data.frame(temp.v=temp.new))
	a.surv.rate <- ifelse(a.surv.rate<0|a.surv.rate>1,0,a.surv.rate)
	return( a.surv.rate ) 
}

## Immature daily emergence rate at different temperature ##
# Correction factor to report the transition time 'acquatic stage' to adult to pupae to adult from Delatte et al. 2009: ((50/100)+(77.5/89.9)+(76.3/93.8)+(67.5/90.9)+(2.5/8.3))/5
i.emer_rate.f <- function(temp.new) {
	i.em_prop.v <- c(0.9,0.91,0.92,0.98,0.98,0.98,0.99,0.98,0.85,0.5,0.5)
	i.em_rate.v <- -log(1-i.em_prop.v)
	i.em_leng.v <- c(61.7,39.7,84.4,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3)
	i.em_drat.v <- i.em_rate.v/i.em_leng.v
	temp.v <- c(14.74,14.84,14.92,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55)
	lm <- lm(i.em_drat.v ~ poly(temp.v,2))
	i.emer.pred <- predict(lm,data.frame(temp.v=temp.new))
	i.emer.pred <- ifelse(i.emer.pred<0,0,i.emer.pred)
	return( i.emer.pred ) 
}
#plot(-15:50,i.emer_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(temp.v,i.em_drat.v)

# png("~/i_emer_rate.png",heigh=480*3,width=480*3,res=300)
# plot(-15:50,i.emer_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,i.emer_rate.f(-15:50),col="blue",ylim=c(0,1))
# legend(x=-10,y=1,c("aegypti","albopictus"),cex=.8,col=c("red","blue"),pch=c(1,2))
# dev.off()

## Immature daily survival rate at different temperature ##
i.surv_rate.f <- function(temp.new) {
	surv.v <- c(6.2,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.7,3.5,0)
	temp.v <- c(10,10.38,10.45,10.45,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64,41)
	lm <- lm(surv.v ~ poly(temp.v,2))
	i.surv.pred <- 1-1/predict(lm,newdata=data.frame(temp.v=temp.new))
	i.surv.pred <- ifelse(i.surv.pred>1|i.surv.pred<0,0,i.surv.pred)
	i.surv.pred[temp.new>45] <- 0
	return( i.surv.pred ) 
}
# plot(seq(-20,50,0.1),i.surv_rate.f(seq(-20,50,0.1)),col="red",type="l",ylim=c(0,1))
# points(1-1/surv.v~temp.v,col="blue")

# png("~/i_surv_rate.png",heigh=480*3,width=480*3,res=300)
# plot(-15:50,i.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,i.surv_rate.f(-15:50),col="blue",ylim=c(0,1))
# legend(x=-10,y=1,c("aegypti","albopictus"),cex=.8,col=c("red","blue"),pch=c(1,2))
# dev.off()

## Immature Density-dependent mortality from Hancock et al. 2009
# Output is an exponential model for daily mortality probability at different larval density
i.dens.v =  c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)*2
i.surv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
i.sur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
	37.12,29.73,40,43,22,51,50,31)
# Transform survival proportion to multidays mortality rate
i.mort_rate.v = -log(i.surv_prop.v)
# Transform multiday survival rate to daily survival rate
i.dmort_rate.v = i.mort_rate.v/i.sur_dur.v
# Fit a model for daily survival rate using density as predictor
i.ddmort_rate.m <- lm(log(i.dmort_rate.v) ~ i.dens.v)

#plot(c(1,100,1000,5000,10000,50000),1-exp(-exp(predict(i.ddmort_rate.m,data.frame(i.dens.v=c(1,100,1000,5000,10000,50000))))),type="l")

## Egg hatching rate
## This rate decides embryonated eggs which hatch or stay from Soares-Pinheiro et al. 2016.
e.betap <- epi.betabuster(mode=0.076, conf=0.95, greaterthan=TRUE, x=0.023)
e.hatch_rate.f <- function(temp.new){
	tmp.p <- rbeta(length(temp.new),e.betap$shape1, e.betap$shape2)
	tmp.p[temp.new<10|temp.new>40] <- 0
	return( tmp.p )
}
#plot(-15:50,e.hatch_rate.f(-15:50),col="red",type="l")

# png("~/e_hatch_rate.png",heigh=480*3,width=480*3,res=300)
# plot(-15:50,e.hatch_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,e.hatch_rate.f(-15:50),col="blue",ylim=c(0,1))
# legend(x=-10,y=1,c("aegypti","albopictus"),cex=.8,col=c("red","blue"),pch=c(1,2))
# dev.off()

## Egg hatching rate: from Thomas et al. 2012 (Fig. 2) and Eisen et al. 2014 (Fig.1)
e.surv_rate.f <- function(temp.new) {
	surv.v <- c(0,0,0,0,0,0,0.4,0.6,0.78,0.81,0.88,0.95,0.96,0.91,0.93,0.83,0.90,0.48)
	temp.v <- c(-17,-15,-12,-10,-7,-5,-2,0,15.6,16,21,22,25.0,26.7,28.0,31.0,32.0,35.0)
	lm <- lm(surv.v ~ poly(temp.v,5))
	e.surv.pred <- predict(lm,data.frame(temp.v=temp.new))
	e.surv.pred <- ifelse(e.surv.pred<0,0,e.surv.pred)
	e.surv.pred[temp.new<= -5] <- 0
	return( e.surv.pred )
}
# plot(-15:50,e.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(surv.v~temp.v,col="blue")

# png("~/e_surv_rate.png",heigh=480*3,width=480*3,res=300)
# plot(-15:50,e.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,e.surv_rate.f(-15:50),col="blue",ylim=c(0,1))
# lines(-15:50,d.surv_rate.f(-15:50),col="green",ylim=c(0,1))
# legend(x=-15,y=1,c("aegypti","albopictus","albopictus dia"),cex=.8,col=c("red","blue","green"),pch=c(1,2,3))
# dev.off()

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)