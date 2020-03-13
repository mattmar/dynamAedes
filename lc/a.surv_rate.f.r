## Adult female daily mortality rate at different temperatures ##
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 4
a.mort_rate.f <- function(dt) {

	survd <- c(13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07)
	st <- c(10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(survd ~ poly(st,4))
	pred_a_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_a_duration<-ifelse(pred_a_duration<1,0.1,pred_a_duration)
	pred_a_rate <- 1/pred_a_duration[,1]
	return( pred_a_rate ) 

}

# su <- a.mort_rate.f(-20:40)
# plot(-20:40,1-(1-exp(-su)))
# lines(st,survd)