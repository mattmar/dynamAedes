##Very basic demographic model for an organism with three life stages (egg, immature [instars and pupa], adult).
#This simulates a single pixel with n introduced individuals. Fixed daily survival and hatching prob.
#s.= surviving, c.=changing stage,.nc= non-changing stage.
#.p = probability, .n=number, .v=vector

#Vector of propagules to initiate the life cycle
i.egg.n <- 20 #n of introduced eggs
i.imm.n <-  0 #n of introduced imamtures
i.adu.n <-  0 #n of introduced adults
days    <-100 #n of day to simulate

#output dataframe where to store numbers of egg, imm and adu per day
res_dt <- as.data.frame(matrix(NA,nrow=days,ncol=3))
names(res_dt) <- c("egg","immature","adult")

for (day in 1:days) {

#Egg compartment
	s.egg.p <- 0.99 #egg survival probability
	c.egg.p <- 0.10 #egg hatching probability
	egg.n   <- ifelse(day==1,i.egg.n,egg.n)
	res_dt[day,1] <- egg.n

#Random draw from binomial distribution to find how many egg die, just survive or survive and hatch
	egg_fate.v <- rmultinom(1, egg.n, c((1-s.egg.p), s.egg.p, c.egg.p))
	nc.egg.n <- egg_fate.v[2,] #surviving eggs which did not hatch
	c.egg.n  <- egg_fate.v[3,] #eggs hatched

#Immature compartment
	s.imm.p <- 0.95 #immature survival probability
	c.imm.p <- 0.10 #immature emerging probability
	nc.imm.n<- ifelse(day==1,i.imm.n,nc.imm.n)
	imm.n   <- c.egg.n+nc.imm.n #updated n of immatures
	res_dt[day,2] <- imm.n
	
#Random draw from binomial distribution to find how many immature die, just survive or survive and emerge
	imm_fate.v <- rmultinom(1, imm.n,c((1-s.imm.p), s.imm.p, c.imm.p))
	nc.imm.n <- imm_fate.v[2,] #immature surviving not emerging
	c.imm.n  <- imm_fate.v[3,] #immature emerging

#Adult compartment
	s.adu.p <- 0.95 #female survival probability
	l.adu.p <- 0.10 #female laying probability
	s.adu.n <- ifelse(day==1,i.adu.n,adu.n)
#Remove males adult from the population using a random binomial draw	
	adu.n <- as.integer(rbinom(1,c.imm.n + s.adu.n,0.5))
	
	res_dt[day,3]  <- adu.n
#Find how many females are ready to lay eggs using a binomial draw	
	l.adu.n <- rbinom(1,adu.n, l.adu.p)

#Random draw from binomial distribution to find how many adults die or survive
	adu_fate.v <- rmultinom(1, adu.n, c((1-s.adu.p), s.adu.p))
	s.adu.n <- adu_fate.v[2,] 

#Find how many of the surviving females lay eggs using a Poisson distribution centred on 40
	l.egg.n <- rpois(l.adu.n,40)
	l.egg.n <-ifelse(length(l.egg.n)==0,0,l.egg.n)

#Updated egg reservoir to start a new "day"
	egg.n <- nc.egg.n+l.egg.n
}

#Plot results
res_lf<-cbind.data.frame(melt(res_dt),days=rep(1:100,3))

ggplot(res_lf, aes(x=days,y=value,col=variable)) + geom_line()
