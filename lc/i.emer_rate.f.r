## Immature daily emergence rate at different temperature ##
## Data taken from https://academic.oup.com/jme/article/46/1/33/902827, Table 10.

# Message to be printed
cat(paste("Imature emergence rate\n"))

#Function
i.emer_rate.f <- function(dt) {

	emerd <- c(1000,61.7,39.7,84.4,15.0,16.0,16.5,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3,1000)
	st <- c(10.45,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,45)
	model <- lm(emerd ~ poly(st,7)) #Sevent polynomial
	pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	pred_emer_rate <- 1/pred_emer_time[,1]
	pred_emer_rate <- ifelse(pred_emer_rate<0,100,pred_emer_rate)
	return( pred_emer_rate ) 

}

# em <- i.emer_rate.f(temp1)
# plot(exp(-(1/emerd))~st,ylim=c(0,1))
# lines(exp(-pred_emer_rate)~dt,col="red")
# plot(pred_emer_time[,1]~dt,ylim=c(0,60))