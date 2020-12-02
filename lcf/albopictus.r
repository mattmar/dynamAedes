#Temperature-dependent life cycle functions for `Aedes albopictus`

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
a.gono_rate.f <- function(temp.new){
  gono.v <- c(35,8.1,4.5,3.5,4.4,10)
  temp.v <- c(10,20,25,30,35,40)
  lm <- lm(gono.v ~ poly(temp.v,4))
  a.gono.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
  a.gono.pred <- ifelse(a.gono.pred<0,0,a.gono.pred)
  return( a.gono.pred )
}
# plot(-15:55,a.gono_rate.f(-15:55),col="red",type="l")
# lines(-15:55,a.gono_rate.f(-15:55),col="blue")
#points(1/gono.v~temp.v,col="blue")

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
a.ovi_rate.f <- function(temp.new){
  ovi.v <-c(50.8, 65.3, 74.2, 48.7)
  temp.v <-c(20,25,30,35)
  lm <- lm(ovi.v ~ poly(temp.v,2))
  a.ovi.pred <- predict(lm,data.frame(temp.v=temp.new))
  a.ovi.pred <- ifelse(a.ovi.pred<0,0,a.ovi.pred)
  return( a.ovi.pred )
}
# plot(-15:50,a.ovi_rate.f(-15:50),col="red",type="l",ylim=c(0,100))
# lines(-15:50,a.ovi_rate.f(-15:50),col="blue")

## Adult female daily mortality rate at different temperatures ##
## Data taken from:
a.surv_rate.f <- function(temp.new){
  a=0.677; b=20.9; c=-13.2
  a.surv.rate <- sapply(temp.new, function(x) {
    if( x>0 ){
      a*exp(-0.5*((x-b)/c)^6)*(x^0.1)
    } else {
      a*exp(-0.5*((x-b)/c)^6)
    }
  })
  return( a.surv.rate )  
}
#plot(-15:50,a.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,a.surv_rate.f(-15:50),col="blue")

#emergence rate immature -> adult: from Tab. 1 in Delatte et al. (2009), column Pupae-adult
i.emer_rate.f <- function(temp.new){
  i.em_prop.v <- c(0,83.3,89.9,93.8,90.0,8.3,0)/100
  i.em_rate.v <- -log(1-i.em_prop.v) 
  i.em_leng.v <- c(20,8.7,4.1,2.7,1.9,1.7,20)
  i.em_drat.v <- i.em_rate.v/i.em_leng.v
  temp.v <- c(-15,15,20,25,30,35,55)
  lm <- lm(i.em_drat.v ~ poly(temp.v,6))
  i.emer.pred <- predict(lm,data.frame(temp.v=temp.new))
  i.emer.pred <- ifelse(i.emer.pred<0,0,i.emer.pred)
  return( 1-exp(-i.emer.pred) )
}
# plot(-15:50,1-exp(-i.emer_rate.f(-15:50)),col="red",type="l",ylim=c(0,1))
# points(1-exp(-i.em_drat.v)~temp.v,col="blue")
# lines(-15:50,i.emer_rate.f(-15:50),col="blue")

## Immature daily survival rate at different temperature ##
## Data taken from:

i.surv_rate.f <- function(dt){
  a=0.977 
  b=20.8 
  c=-12.6
  i.surv.pred <- a*exp(-0.5*((dt-b)/c)^6)
  return( i.surv.pred ) 
}
# plot(-15:50,i.surv_rate.f(-15:50),col="red",type="l")
# lines(-15:50,i.surv_rate.f(-15:50),col="blue")

## Density-dependent immature mortality, fit data from Hancock et al. 2009
# Output is an exponential model for daily mortality probability at different larval density
i.dens.v =  c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)*2
i.surv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
i.sur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
  37.12,29.73,40,43,22,51,50,31)
# Transform survival proportion to multidays mortality rate
i.mort_rate.v = -log(i.surv_prop.v)
# Transform multiday survival rate to daily survival rate
i.dmort_rate.v = i.mort_rate.v/i.sur_dur.v
# Fit a model for daily survival rate using density as predictor
i.ddmort_rate.m <- lm(log(i.dmort_rate.v) ~ i.dens.v)

#plot(c(1,100,1000,5000,10000,50000),1-exp(-exp(predict(i.ddmort_rate.m,data.frame(i.dens.v=c(1,100,1000,5000,10000,50000))))),type="l")

## Egg hatching rate: from Tab. 1 in Delatte et al. (2009), column Egg-L1
## This rate decides embryonated eggs which hatch or stay
e.hatch_rate.f <- function(temp.new){
  e.hatch_leng.v <- c(100,7.4,2.9,4.5,6.7,7.1,100) 
  temp.v <- c(0,15,20,25,30,35,55)
  lm <- lm(e.hatch_leng.v ~ poly(temp.v,4))
  e.hatch.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
  e.hatch.pred[temp.new<10] <- 0
  return( e.hatch.pred )
}
# plot(-15:50,e.hatch_rate.f(-15:50),col="red",type="l",ylim=c(0,0.50))
# points(1/e.hatch_leng.v~temp.v)
# points(-15:50,e.hatch_rate.f(-15:50),col="blue")

#From Metelmann et al. (2019)
## These probabilities decide which eggs die or survive
e.surv_rate.f <- function(dt){
  a=0.95 
  b=18.8 
  c=-21.53
  e.surv.pred=a*exp(-0.5*((dt-b)/c)^6)
  return( e.surv.pred ) 
}
#plot(-10:50,e.surv_rate.f(-10:50),col="red",ylim=c(0,1),type="l")

d.surv_rate.f <- function(dt){
  ed_surv_bl=1
  a=0.99
  b=11.68
  c=-15.67
  d.surv.pred= ed_surv_bl*a*exp(-0.5*((dt-b)/c)^6)
  return( d.surv.pred )
}
#plot(-20:50,d.surv_rate.f(-20:50),col="red",ylim=c(0,1),type="l")
#d.surv_rate.f(-10)
# points(-10:50,e.surv_rate.f(-10:50),col="blue")

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)