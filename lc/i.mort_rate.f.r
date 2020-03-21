## Immature daily survival rate at different temperature ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 4.
i.mort_rate.f <- function(dt) {

	survd <- c(0,6.2,10.3,12.2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.7,3.5,0)
	st <- c(0,10,10,10,10,10.38,10.45,10.45,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64,50)
	model <- lm(survd ~ poly(st,3))
	pred_i_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_i_duration[,1]<-ifelse(pred_i_duration[,1]<0,000.1,pred_i_duration[,1])
	pred_i_rate <- 1/pred_i_duration[,1]
	return( pred_i_rate ) 

}

# su <- i.mort_rate.f(-10:45)
# plot(-10:45,1-(1-(exp(-su))),type="l")
# lines(st,1-(1/survd),col="red")