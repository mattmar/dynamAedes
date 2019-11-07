## Daily rate of gonotrophic cycle, i.e., rate of adult females completing the cycle from bloodmeal to the first day of oviposition.
## Data taken from https://academic.oup.com/jme/article/46/1/33/902827, Table 5.

# Message to be printed
cat(paste("Gonotrophic cycle T-dependent rate\n"))

# Function
a.gono_rate.f <- function(dt) {

	days <- c(60,30,15,7,5,5,6,5,6,5,6,10)
	st <- c(10,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,40)
	model <- lm(days ~ poly(st,3)) #Third polynomial function
	pred_gono_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence', level=0.95)
	return(1/pred_gono_duration[,1]) # Return a rate

}

#gp <- a.gono_duration.f(temp1)
#plot(temp1,a.batch.n)