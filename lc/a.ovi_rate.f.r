## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
## Data taken from https://academic.oup.com/jme/article/46/1/33/902827, Table 5.

# Message to be printed
cat(paste("Oviposition rate per female per day\n"))

# Function
a.ovi_rate.f <- function(dt) {

	ovir <- c(0.3548,1.1208,3.3668,3.5931,6.9847,7.5997,9.5762,7.2770,11.224,7.2745)
	st <- c(15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41)
	model <- lm(ovir ~ poly(st,4)) #Forth polynomial
	pred_ovi_rate <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	return(pred_ovi_rate[,1])

}

# ov <- a.ovi_rate.f(temp1)
# plot(temp1,a.batch.n)
