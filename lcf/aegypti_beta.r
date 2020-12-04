#Temperature-dependent life cycle functions for `Aedes aegypti`

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
a.gono_rate.f <- function(temp.new) {
	gono_d <- 1/c(30,15,7,5,4,4,5,5,6,7,35,100)
	temp_d <- c(15,17,20,22,26,28,31,32,33,33,40, 45)
	model <- drm(gono_d ~ temp_d, fct = DRC.beta())
	a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
	return( a.gono.pred ) 
}
#plot(-15:55,a.gono_rate.f(-15:55),col="red",type="l")
# points(1/gono.v~temp.v,col="blue")

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
a.ovi_rate.f <- function(temp.new) {
	ovi_n <- c(0,0,0,0.35,1.12,3.37,3.6,6.98,7.6,9.58,7.28,11.22,7.27,5,0) * 10
	temp_n <- c(0,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,34,45)
	model <- drm(ovi_n ~ temp_n, fct = DRC.beta())
	a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
	return( a.ovi.pred )
}
#plot(-15:55,a.ovi_rate.f(-15:55),col="red",type="l")
# points(ovi_n~temp_n,col="blue")

## Adult female daily mortality rate at different temperatures ##
a.surv_rate.f <- function(temp.new) {
	surv_d <- c(3,13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07,5)
	temp_d <- c(8,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,40)
	model <- drm(surv_d ~ temp_d, fct = DRC.beta())
	a.surv.rate <- predict(model,data.frame(temp.v=temp.new))
	a.surv.rate <- 1-1/ifelse(a.surv.rate<1,1,a.surv.rate)
	return( a.surv.rate ) 
}
#plot(-15:55,a.surv_rate.f(-15:55),col="red",type="l")
# points(1-1/surv_d~temp_d,col="blue")
# cbind(-15:55,round(a.surv_rate.f(-15:55),2))

## Immature daily emergence rate at different temperature ##
i.emer_rate.f <- function(temp.new) {
	eme_d <- 1/(c(500,61.7,39.7,84.4,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3,30)*0.6439)
	temp_d <- c(12,14.74,14.84,14.92,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,45)
	model <- drm(eme_d ~ temp_d, fct = DRC.beta())
	i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
	return( i.emer.pred ) 
}
#plot(-15:50,i.emer_rate.f(-15:50),col="red",type="l")
# points(temp_d,eme_d)
# cbind(-15:55,round(i.emer_rate.f(-15:55),2))

## Immature daily survival rate at different temperature ##
i.surv_rate.f <- function(temp.new) {
	surv_d <- c(6.2,10.3,12,2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.7,3.5,1)
	temp_d <- c(10,10,10,10,10.38,10.45,10.45,10,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64,45)
	model <- drm(surv_d ~ temp_d, fct = DRC.beta())
	i.surv.pred <- predict(model,data.frame(temp.v=temp.new))
	i.surv.pred <- 1-1/ifelse(i.surv.pred<1,1,i.surv.pred)
	return( i.surv.pred ) 
}
# plot(-15:50,i.surv_rate.f(-15:50),col="red",type="l")
# points(1-1/surv.v~temp.v,col="blue")
#cbind(-15:55,round(i.surv_rate.f(-15:55),2))

## Immature Density-dependent mortality from Hancock et al. 2009
# Output is an exponential model for daily mortality probability at different larval densities; the density refers to 2L of water (in a 5L container)
i.dens.v = c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)
i.surv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
i.sur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
	37.12,29.73,40,43,22,51,50,31)
# Transform survival proportion to multidays mortality rate
i.mort_rate.v = -log(i.surv_prop.v)
# Transform multiday mortality rate to daily mortality rate
i.dmort_rate.v = i.mort_rate.v/i.sur_dur.v
# Fit a model for daily mortality rate using just density as a predictor
i.ddmort_rate.m <- lm(log(i.dmort_rate.v) ~ i.dens.v)

#plot(c(1,100,1000,5000,10000,50000),1-exp(-exp(predict(i.ddmort_rate.m,data.frame(i.dens.v=c(1,100,1000,5000,10000,50000))))),type="l")

## Egg hatching rate
## This rate decides embryonated eggs which hatch or stay
e.hatch_rate.f <- function(temp.new){
  hatc_d <- c(0,  0.025*2,1,0.98, 0.73, 0.30, 0.037, 0.016)/2 
  temp_d <- c(7, 12,   19,  24.5, 26.5, 29.5, 32.5,    34.5)
  model <- drm(hatc_d ~ temp_d, fct = DRC.beta())
  e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
  return( e.hatch.pred )
}
# plot(-15:50,e.hatch_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(hatc_d~temp_d,col="blue")
# cbind(-15:55,round(e.hatch_rate.f(-15:55),2))

## Egg hatching rate: from Thomas et al. 2012 (Fig. 2) and Eisen et al. 2014 (Fig.1)
e.surv_rate.f <- function(temp.new) {
	surv_r <- c(0,0,0,0,0,0,0.4,0.6,0.78,0.81,0.88,0.95,0.96,0.91,0.93,0.83,0.90,0.48,0)
	temp_r <- c(-17,-15,-12,-10,-7,-5,-2,0,15.6,16,21,22,25.0,26.7,28.0,31.0,32.0,35.0,45.0)
	model <- drm(surv_r ~ temp_r, fct = DRC.beta())
	e.surv.pred <- predict(model,data.frame(temp.v=temp.new))
	return( e.surv.pred )
}
# plot(-15:50,e.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(surv.v~temp.v,col="blue")
# cbind(-15:55,round(e.surv_rate.f(-15:55),2))

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)