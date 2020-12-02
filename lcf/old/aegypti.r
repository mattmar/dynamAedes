#Temperature-dependent life cycle functions for `Aedes aegypti`

## Immature daily emergence rate at different temperature ##
## Data taken from 
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 10
i.emer_rate.f <- function(dt) {
	emerd <- c(61.7,39.7,84.4,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3)
	st <- c(14.74,14.84,14.92,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55)
	model <- lm(emerd ~ poly(st,2))
	pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
	pred_emer_time <- 1-ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])/100
	return( 1-pred_emer_time ) 
}

## Immature daily survival rate at different temperature ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 4.
i.mort_rate.f <- function(dt) {
	survd <- c(0,6.2,10.3,12.2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.7,3.5,0)
	st <- c(0,10,10,10,10,10.38,10.45,10.45,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64,50)
	model <- lm(survd ~ poly(st,3))
	pred_i_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_i_duration[,1] <- ifelse(pred_i_duration[,1]<0,000.1,pred_i_duration[,1])
	pred_i_rate <- 1/pred_i_duration[,1]
	return( pred_i_rate ) 
}

## Adult female daily mortality rate at different temperatures ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 4
a.mort_rate.f <- function(dt) {
	survd <- c(0,0,13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07)
	st <- c(-20,5,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(survd ~ poly(st,3))
	pred_a_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_a_duration <- ifelse(pred_a_duration[,1]<1,0,pred_a_duration[,1])
	pred_a_rate <- 1/pred_a_duration
	pred_a_rate[pred_a_rate==Inf] <- 1
	return( pred_a_rate ) 
}

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
# This process has a minimum duration of 5 days, which consiste of 24h to mature gonads, 3 days for gonadotrophic cycle.
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 5
a.gono_rate.f <- function(dt) {
	gpd <- c(30,15,7,5,5,6,5,6,5,6)
	st <- c(15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(gpd ~ poly(st,6))
	pred_gono_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
	pred_gono_duration[,1] <- ifelse(pred_gono_duration[,1]<0,0,pred_gono_duration[,1])
	pred_gono_rate <- 1/pred_gono_duration[,1]
	return( pred_gono_rate ) 
}

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 5
a.ovi_rate.f <- function(dt) {
	ovir <- c(0,0,0.3548,1.1208,3.3668,3.5931,6.9847,7.5997,9.5762,7.2770,11.224,7.2745)
	ovir <- ovir * 10
	st <- c(-10,10.54,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(ovir ~ poly(st,3))
	pred_ovi_rate <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_ovi_rate <- ifelse(pred_ovi_rate[,1]<0,0,pred_ovi_rate[,1])
	return( pred_ovi_rate )
}

## Uniform rate of egg survival; set to have a uniform p=0.95 as it can be assumed that egg survival is independent from temperature.
# e.surv_rate.f <- function(dt){
# 	return( -log(1-0.98) )
# }
e.mort_rate.f <- function(dt) {
	survr <- c(0,0,0,0,0,0.4,0.6,0.99,0.99,1.00,1.00,0.70,0)
	st <- c(-15,-12,-10,-7,-5,-2,0,10,20,25,30,40,50)
	model <- lm(survr ~ poly(st,6))
	pred_surv_p <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_surv_p <- ifelse(pred_surv_p[,1]<0,0,pred_surv_p[,1])
	pred_surv_p <- ifelse(pred_surv_p>=0.99,0.99,pred_surv_p)
	return( 1-pred_surv_p )
}
# plot(-15:55,e.surv_rate.f(-15:55),col="red")
# points(-15:55,1-exp(-e.surv_rate.f(-15:55)),col="blue")

## Egg hatching rate (from DOI: 10.1590/S0037-86822012000200007)
e.betap <- epi.betabuster(mode=0.076, conf=0.95, greaterthan=TRUE, x=0.023)
e.hatch_rate.f <- function(dt){
	tmp.p <- rbeta(length(dt),e.betap$shape1, e.betap$shape2)
	return( 1-tmp.p )
}

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)