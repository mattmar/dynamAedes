## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 5
a.ovi_rate.f <- function(dt) {

	ovir <- c(0.3548,1.1208,3.3668,3.5931,6.9847,7.5997,9.5762,7.2770,11.224,7.2745)
	st <- c(15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(ovir ~ poly(st,2))
	pred_ovi_rate <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_ovi_rate<-ifelse(pred_ovi_rate<0,0,pred_ovi_rate)
	return(pred_ovi_rate[,1])

}

# ovr <- a.ovi_rate.f(1:50)
# plot(1:50,ovr)
# lines(st,ovir,col="red")