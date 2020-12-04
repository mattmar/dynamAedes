#Temperature-dependent life cycle functions for `Aedes albopictus`

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
a.gono_rate.f <- function(temp.new){
  gono_d <- 1/c(8.1, 4.5, 3.5, 4.4,100)
  temp_d <-   c(20,  25,  30,  35, 45)
  model <- drm(gono_d ~ temp_d, fct = DRC.beta())
  a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
  return( a.gono.pred )
}
# plot(-10:55,a.gono_rate.f(-10:55),col="red",type="l")
# lines(-15:55,a.gono_rate.f(-15:55),col="blue")
# points(1/gono.v~temp.v,col="blue")
# cbind(-15:55,round(a.gono_rate.f(-15:55),2))

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
a.ovi_rate.f <- function(temp.new){
  ovi_n <-  c(0, 50.8, 65.3, 74.2, 48.7, 0)
  temp_n <- c(5, 20,   25,   30,   35, 45)
  model <- drm(ovi_n ~ temp_n, fct = DRC.beta())
  a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
  return( a.ovi.pred )
}
#plot(-15:50,a.ovi_rate.f(-15:50),col="red",type="l",ylim=c(0,100))
# lines(-15:50,a.ovi_rate.f(-15:50),col="blue")
#cbind(-15:55,round(a.ovi_rate.f(-15:55),2))

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
# plot(-15:50,a.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# lines(-15:50,a.surv_rate.f(-15:50),col="blue")
# cbind(-15:55,round(a.surv_rate.f(-15:55),2))

#emergence rate immature -> adult: from Tab. 1 in Delatte et al. (2009), column Pupae-adult
i.emer_rate.f <- function(temp.new){
  eme_d <- 1/c(100,8.7,4.1,2.7,1.9,1.7,5)
  temp_d <- c(5,15,20,25,30,35,40)
  model <- drm(eme_d ~ temp_d, fct = DRC.beta())
  i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
  return( i.emer.pred )
}
#plot(-15:50,1-exp(-i.emer_rate.f(-15:50)),col="red",type="l",ylim=c(0,1))
#points(eme_d~temp_d,col="blue")
# lines(-15:50,i.emer_rate.f(-15:50),col="blue")
# cbind(-15:55,round(i.emer_rate.f(-15:55),2))

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
# cbind(-15:55,round(i.surv_rate.f(-15:55),2))

## Immature Density-dependent mortality from Hancock et al. 2009
# Output is an exponential model for daily mortality probability at different larval densities; the density refers to 2L of water (in a 5L container)
i.dens.v = c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)
i.surv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
i.sur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
  37.12,29.73,40,43,22,51,50,31)
# Transform survival proportion to multidays mortality rate
i.mort_rate.v = -log(i.surv_prop.v)
# Transform multiday mortality rate to daily mortality rate
i.dmort_rate.v = i.mort_rate.v/i.sur_dur.v
# Fit a model for daily mortality rate using just density as a predictor
i.ddmort_rate.m <- lm(log(i.dmort_rate.v) ~ i.dens.v)

#plot(c(1,100,1000,5000,10000,50000),1-exp(-exp(predict(i.ddmort_rate.m,data.frame(i.dens.v=c(1,100,1000,5000,10000,50000))))),type="l")

## Egg hatching rate: from Tab. 1 in Delatte et al. (2009), column Egg-L1
## This rate decides embryonated eggs which hatch or stay
e.hatch_rate.f <- function(temp.new){
  hatc_d <- 1/c(100,11,7.4,2.9,4.5,6.7,7.1,500) 
  temp_d <- c(0,10,15,20,25,30,35,35)
  model <- drm(hatc_d ~ temp_d, fct = DRC.beta())
  e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
  return( e.hatch.pred )
}
# plot(-15:50,e.hatch_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(hatc_d~temp_d)
# points(-15:50,e.hatch_rate.f(-15:50),col="blue")
# cbind(-15:55,round(e.hatch_rate.f(-15:55),2))

#From Metelmann et al. (2019)
## These probabilities decide which eggs die or survive
e.surv_rate.f <- function(dt){
  a=0.955 
  b=18.8 
  c=-21.53
  e.surv.pred=a*exp(-0.5*((dt-b)/c)^6)
  return( e.surv.pred ) 
}
#plot(-10:50,e.surv_rate.f(-10:50),col="red",ylim=c(0,1),type="l")

d.surv_rate.f <- function(dt){
  ed_surv_bl=1
  a=0.955
  b=11.68
  c=-15.67
  d.surv.pred= ed_surv_bl*a*exp(-0.5*((dt-b)/c)^6)
  return( d.surv.pred )
}
# plot(-20:50,d.surv_rate.f(-20:50),col="red",ylim=c(0,1),type="l")
# d.surv_rate.f(-10)
# points(-10:50,e.surv_rate.f(-10:50),col="blue")
# cbind(-15:55,round(d.surv_rate.f(-15:55),2))

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)
#cbind(seq(0,600,10),round(f.adis.p,5))
