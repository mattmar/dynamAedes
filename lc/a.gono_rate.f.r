## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition.
# This process has a minimum duration of 5 days, which consiste of 24h to mature gonads, 3 days for gonadotrophic cycle.
## Data taken from:
## https://www.cambridge.org/core/journals/epidemiology-and-infection/article/assessing-the-effects-of-temperature-on-the-population-of-aedes-aegypti-the-vector-of-dengue/E2FE126FB84D0DE97A94E68343B4649C/core-reader
## Table 5
a.gono_rate.f <- function(dt) {

	gpd <- c(30,15,7,5,5,6,5,6,5,6)
	st <- c(15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(gpd ~ poly(st,6))
	pred_gono_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
	pred_gono_duration[,1]<-ifelse(pred_gono_duration[,1]<0,0,pred_gono_duration[,1])
	pred_gono_rate <- 1/pred_gono_duration[,1]
	return( pred_gono_rate ) 

}

 # gp <- a.gono_rate.f(-10:45)
 # plot(-10:45,1-exp(-gp),type="l")
 # lines(st,(1/gpd),col="red")