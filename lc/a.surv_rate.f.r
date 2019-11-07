## Adult female daily survival rate at different temperature ##
## Data taken from https://academic.oup.com/jme/article/46/1/33/902827, Table 4.

# Message to be printed
cat(paste("Adult female survival rate\n"))

#Function
a.surv_rate.f <- function(dt) {

	survd <- c(13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,20.07)
	st <- c(10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(survd ~ poly(st,4)) #Forth polynomial
	pred_a_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	return((1/pred_a_duration[,1]))

}

#su <- a.surv_rate.f(temp1)
#plot(temp1,su)
