## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
.a.gono_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		gono_d <- 1/c(30,15,7,5,4,4,5,5,6,7,35,100)
		temp_d <- c(15,17,20,22,26,28,31,32,33,33,40, 45)
		model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
		a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="albopictus") {
		gono_d <- 1/c(8.1, 4.5, 3.5, 4.4,100)
		temp_d <-   c(20,  25,  30,  35, 45)
		model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
		a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="koreicus") {
		gono_d <- 1/c(14.75, 11.5, 9.21, 10.81, 100)
		temp_d <- c(  18,    23,   23,    28,    33)
		model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
		a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
	}else(stop("Species not supported."))
	return(a.gono.pred)
}

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
.a.ovi_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		ovi_n <- c(0,0,0,0.35,1.12,3.37,3.6,6.98,7.6,9.58,7.28,11.22,7.27,5,0) * 10
		temp_n <- c(0,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,34,45)
		model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
		a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="albopictus") {
		ovi_n <-  c(0, 50.8, 65.3, 74.2, 48.7, 0)
		temp_n <- c(5, 20,   25,   30,   35, 45)
		model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
		a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="koreicus") {
		ovi_n <- c(0,50.8, 65.3, 74.2, 48.7,0)/0.7347
		temp_n <-c(0,20, 25,30,35,45)
		model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
		a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
	}else(stop("Species not supported."))
}

## Adult female daily mortality rate at different temperatures ##
.a.surv_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		surv_d <- c(3,13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07,5)
		temp_d <- c(8,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,40)
		model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
		a.surv.rate <- predict(model,data.frame(temp.v=temp.new))
		a.surv.rate <- 1-1/ifelse(a.surv.rate<1,1,a.surv.rate)
	}else if(sp=="albopictus") {
		.a.surv_rate.f <- function(temp.new){
			a=0.677; b=20.9; c=-13.2
			a.surv.rate <- sapply(temp.new, function(x) {
				if( x>0 ){
					a*exp(-0.5*((x-b)/c)^6)*(x^0.1)
				} else {
					a*exp(-0.5*((x-b)/c)^6)
				}
			})
		}
	}else if(sp=="koreicus") {
		surv_d <- c(2,52.33,76.77,66.33,5.87,10)
		temp_d <- c(5,18,23,28,33,35)
		model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
		a.surv.rate <- predict(model,data.frame(temp.v=temp.new))
		a.surv.rate <- 1-1/ifelse(a.surv.rate<1,1,a.surv.rate)
	}else(stop("Species not supported."))
}

## Log-Normal probability density of short active dispersal from Marcantonio et al, 2019; Marini et al. 2019. 0 to 600 m with resolution of 10 m.
.a.a_disp.f <- function(sp, max.a.disp, disp.bins) {
	if(sp=="aegypti"){
		dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.95, sdlog=0.66)
	}else if(sp=="albopictus") {
		dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.54, sdlog=0.58)
	}else if(sp=="koreicus") {
		dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.54, sdlog=0.58)
	}else(stop("Species not supported."))
}

## Immature daily emergence rate at different temperature ##
.i.emer_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		eme_d <- 1/(c(500,61.7,39.7,84.4,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3,30)*0.6439)
		temp_d <- c(12,14.74,14.84,14.92,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,45)
		model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
		i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="albopictus") {
		eme_d <- 1/c(100,8.7,4.1,2.7,1.9,1.7,5)
		temp_d <- c(5,15,20,25,30,35,40)
		model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
		i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="koreicus") {
		eme_d <- 1/c(35, 10.31, 4.19, 3.43,  3.00, 2.07, 1.82, 5, 25)
		temp_d <-  c(8,  13,    18,   23,    23,   28,   33,  35, 40)
		model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
		i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
	}else(stop("Species not supported."))
}

## Immature daily survival rate at different temperature ##
.i.surv_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		surv_d <- c(6.2,10.3,12,2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.7,3.5,1)
		temp_d <- c(10,10,10,10,10.38,10.45,10.45,10,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64,45)
		model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
		i.surv.pred <- predict(model,data.frame(temp.v=temp.new))
		i.surv.pred <- 1-1/ifelse(i.surv.pred<1,1,i.surv.pred)
	}else if(sp=="albopictus") {
		a=0.977
		b=20.8 
		c=-12.6
		i.surv.pred <- a*exp(-0.5*((temp.new-b)/c)^6)
	}else if(sp=="koreicus") {
		surv_d <- c(0, 80, 97.8, 93.5, 95.6, 92.8, 84.7, 0)/100 
		temp_d <- c(4, 8, 13,   18,   23,   28,   33,   45)
		model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
		i.surv.pred <- predict(model,data.frame(temp.v=temp.new))
	}else(stop("Species not supported."))
}

## Immature Density-dependent mortality from Hancock et al. 2009
.i.ddmort_rate.f <- function(dens.new) {
	i.dens.v = c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)
	i.dsurv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
	i.dsur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
		37.12,29.73,40,43,22,51,50,31)
	# Transform survival proportion to multidays mortality rate
	i.dmort_rate.v = -log(i.dsurv_prop.v)
	# Transform multiday mortality rate to daily mortality rate
	i.dmort_rate.v = i.dmort_rate.v/i.dsur_dur.v
	# Fit a model for daily mortality rate using just density as a predictor
	i.dmort.lm <- lm(log(i.dmort_rate.v) ~ i.dens.v)
	predict(i.dmort.lm, dens.new)
}

## Egg hatching rate
## This rate decides embryonated eggs which hatch or stay
.e.hatch_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		hatc_d <- c(0,  0.025*2,1,0.98, 0.73, 0.30, 0.037, 0.016)/2 
		temp_d <- c(7, 12,   19,  24.5, 26.5, 29.5, 32.5,    34.5)
		model <- drm(hatc_d ~ temp_d, fct = .DRC.beta())
		e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="albopictus") {
		hatc_d <- 1/c(100,11,7.4,2.9,4.5,6.7,7.1,500) 
		temp_d <- c(0,10,15,20,25,30,35,35)
		model <- drm(hatc_d ~ temp_d, fct = .DRC.beta())
		e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="koreicus") {
		hatc_d <- 1/c(100,2.45,1.35,1.07,1.08,2.04,100)/2
		temp_d <- c(7,8,13,23,28,33,35)
		model <- drm(hatc_d ~ temp_d, fct = .DRC.beta())
		e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
	}else(stop("Species not supported."))
}

## Egg hatching rate: from Thomas et al. 2012 (Fig. 2) and Eisen et al. 2014 (Fig.1)
.e.surv_rate.f <- function(temp.new, sp) {
	if(sp=="aegypti") {
		surv_r <- c(0,0,0,0,0,0,0.4,0.6,0.78,0.81,0.88,0.95,0.96,0.91,0.93,0.83,0.90,0.48,0)
		temp_r <- c(-17,-15,-12,-10,-7,-5,-2,0,15.6,16,21,22,25.0,26.7,28.0,31.0,32.0,35.0,45.0)
		model <- drm(surv_r ~ temp_r, fct = .DRC.beta())
		e.surv.pred <- predict(model,data.frame(temp.v=temp.new))
	}else if(sp=="albopictus") {
		a=0.955 
		b=18.8 
		c=-21.53
		e.surv.pred=a*exp(-0.5*((temp.new-b)/c)^6)
	}else if(sp=="koreicus") {
		ed_surv_bl=1
		a=0.98
		b=15.68
		c=-17.67
		d.surv.pred= ed_surv_bl*a*exp(-0.5*((temp.new-b)/c)^6)
	}else(stop("Species not supported."))
}

## Diapausing hatching rate: from Thomas et al. 2012 (Fig. 2) and Eisen et al. 2014 (Fig.1)
.d.surv_rate.f <- function(temp.new){
	ed_surv_bl=1
	a=0.955
	b=11.68
	c=-15.67
	d.surv.pred=ed_surv_bl*a*exp(-0.5*((temp.new-b)/c)^6)
	return( d.surv.pred )
}