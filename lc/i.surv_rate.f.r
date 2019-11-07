## Immature daily survival rate at different temperature ##
## Data taken from https://academic.oup.com/jme/article/46/1/33/902827, Table 4.

# Message to be printed
cat(paste("Immature survival rate\n"))

#Function
i.surv_rate.f <- function(dt) {

	survd <- c(6.2,10.3,12.2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.8,15.8,15.7,9.5,16.8,8.6,14.2,10.0,3.1,3.8,3.5)
	st <- c(10,10,10,10,10.38,10.45,10.45,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,39.95,40.16,40.64)
	model <- lm(survd ~ poly(st,4)) #Forth polynomial
	pred_i_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	return(1/pred_i_duration[,1])

}

# su <- i.surv_rate.f(10:40)
# plot(10:40,exp(-su),ylim=c(0,1),xlim=c(10,40),type="l")
# points(st,exp(-1/survd),col="red")
