## Immature daily emergence rate at different temperature ##
## Data taken from 
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 10
i.emer_rate.f <- function(dt) {

	emerd <- c(61.7,39.7,84.4,15.0,16.0,16.5,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3)
	st <- c(14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55)
	model <- lm(emerd ~ poly(st,2))
	pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
	pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,10e6,pred_emer_time[,1])
	pred_emer_rate <- 1/pred_emer_time[,1]
	return( pred_emer_rate ) 

}

# em <- i.emer_rate.f(-20:50)
# plot(-20:50,1-exp(-em))
# lines(st,1/emerd,col="red")