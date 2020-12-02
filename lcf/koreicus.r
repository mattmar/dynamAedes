#Temperature-dependent life cycle functions for `Aedes koreicus`

# Gonotrophic cycle. Marini et al. (2019) and Ciocchetta et al. (2017).
a.gono_rate.f <- function(temp.new){
  gono.v <-c(100, 14.75, 11.5, 9.21, 10.81, 30, 100)
  temp.v <-c(0,   18,    23,  23,   28,    33, 40)
  lm <- lm(gono.v ~ poly(temp.v,3))
  a.gono.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
  return( a.gono.pred ) 
}
# plot(-10:40,a.gono_rate.f(-10:40),col="red",type="l")
# points(1/gono.v~temp.v,col="blue")

## Oviposition rate, i.e., number of eggs laid per female/day ##
a.ovi_rate.f <- function(temp.new){
  ovi.v <- c(0,50.8, 65.3, 74.2, 48.7,0)/0.653
  temp.v <-c(0,20, 25,30,35,40)
  lm <- lm(ovi.v ~ poly(temp.v,4))
  a.ovi.pred <- predict(lm,data.frame(temp.v=temp.new))
  a.ovi.pred <- ifelse(a.ovi.pred<0,0,a.ovi.pred)
  a.ovi.pred[temp.new<0] <- 0
  return( a.ovi.pred )
}
# plot(-15:50,a.ovi_rate.f(-15:50),col="red",type="l",ylim=c(0,120))
# points(ovi.v~temp.v,col="blue")

# a.ovi_rate.f <- function(temp.new){
#   a.ovi.pred <- rnorm(length(temp.new),mean=100.4,sd=35.95)
#   return( a.ovi.pred )
# }

#From: Marini et al. (2019)
a.surv_rate.f <- function(temp.new) {
  surv.time <- c(1,2,52.33,46.77,66.33,5.87,2)
  temp.v <- c(-5,5,18,23,28,33,39)
  lm <- lm(surv.time ~ poly(temp.v,3))
  a.surv.rate <- 1-1/predict(lm,data.frame(temp.v=temp.new))
  a.surv.rate[temp.new<= -0|temp.new>38] <-0
  return( a.surv.rate ) 
}
# plot(-15:50,a.surv_rate.f(-15:50),col="red",type="l",ylim=c(0,1))
# points(temp.v,1-1/surv.time,col="blue")

#emergence rate immature
#development: duration of development (number of days) from Tab. 2
#Immature: average for colum Pupae
# i.emer_rate.f=function(temp.new){
#   i.em_leng.v <- c(95,25, 9.062, 4.542, 2.95, 2.42,2.486,20 )
#   temp.v <- c(0,8,13,18,23,28,33,40)
#   lm <- lm(i.em_leng.v ~ poly(temp.v,7))
#   i.emer.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
#   i.emer.pred[temp.new<0] <- 0
#   return( i.emer.pred )
# }
i.emer_rate.f=function(temp.new){
  i.em_leng.v <- c(25, 10.31, 4.19, 3.43,  3.00, 2.07, 1.82, 5)
  temp.v <-      c(8,  13,    18,   23,    23,   28,   33,   40)
  lm <- lm(i.em_leng.v ~ poly(temp.v,4))
  i.emer.pred <- 1/predict(lm,data.frame(temp.v=temp.new))
  i.emer.pred[temp.new<0] <- 0
  return( i.emer.pred )
}
# plot(-15:50,i.emer_rate.f(-15:50),col="red",type="l")
# points(1/i.em_leng.v~temp.v,col="blue")

#immature survival rate. Marini et al. (2019)
i.surv_rate.f <- function(temp.new){
  surv.v <-c(0 , 97.8, 93.5, 95.6, 92.8, 84.7, 0)/100 
  temp.v <-c(4,  13,   18,   23,   28,   33,   40)
  lm <- lm(surv.v ~ poly(temp.v,4))
  i.surv.pred <- predict(lm,data.frame(temp.v=temp.new))
  i.surv.pred <- ifelse(i.surv.pred<0,0,i.surv.pred)
  return( i.surv.pred ) 
}
# plot(-15:50,i.surv_rate.f(-15:50),col="red",type="l")
# points(surv.v~temp.v,col="blue")

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

#eggs hatching rate
e.hatch_rate.f <- function(temp.new){
  emer.l <- c(100,2.45,1.35,1.07,1,1.04,100)
  temp.v <- c(-10,8,13,23,28,33,60)
  lm <- lm(emer.l ~ poly(temp.v,6))
  e.hatch_rate.v <- 1/predict(lm,data.frame(temp.v=temp.new))
  e.hatch_rate.v[temp.new<10] <- 0
  e.hatch_rate.v[e.hatch_rate.v>1] <- 1
  return( e.hatch_rate.v ) 
}
#plot(-10:50,e.hatch_rate.f(-10:50),col="red",type="l")
# points(1/emer.l~temp.v,col="blue")

#eggs survival rate.
e.surv_rate.f <- function(temp.new){
  surv.v <- c(0, 20,  50.50, 53.75, 51.00, 57.25,0)/100
  temp.v <- c(-12, 0,  13,    23,    28,    33,40)
  lm <- lm(surv.v ~ poly(temp.v,4))
  e.surv.pred <- predict(lm,data.frame(temp.v=temp.new))
  e.surv.pred <- ifelse(e.surv.pred<0,0,e.surv.pred)
  return( e.surv.pred )
}
# plot(-10:50,e.surv_rate.f(-10:50),col="red",type="l")
# points(surv.v~temp.v,col="blue")

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)