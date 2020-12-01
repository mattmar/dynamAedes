#Temperature-dependent life cycle functions for `Aedes koreicus`

# Gonotrophic cycle. Marini et al. (2019) and Ciocchetta et al. (2017).
a.gono_rate.f <- function(temp.new){
  gono_d <- 1/c(14.75, 11.5, 9.21, 10.81, 100)
  temp_d <- c(  18,    23,   23,    28,    33)
  model <- drm(gono_d ~ temp_d, fct = DRC.beta())
  a.gono.pred <- predict(model,data.frame(temp.v=temp.new))
  return( a.gono.pred ) 
}

# plot(-10:40,a.gono_rate.f(-10:40),col="red",type="l")
# points(gono_d~temp_d,col="blue")
# cbind(-15:55,round(a.gono_rate.f(-15:55),2))

## Oviposition rate, i.e., number of eggs laid per female/day ##
a.ovi_rate.f <- function(temp.new){
  ovi_n <- c(0,50.8, 65.3, 74.2, 48.7,0)/0.7347
  temp_n <-c(0,20, 25,30,35,45)
  model <- drm(ovi_n ~ temp_n, fct = DRC.beta())
  a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
  return( a.ovi.pred )
}
# plot(-15:50,a.ovi_rate.f(-15:50),col="red",type="l",ylim=c(0,120))
# points(ovi.v~temp.v,col="blue")
# cbind(-15:55,round(a.ovi_rate.f(-15:55),2))

#From: Marini et al. (2019)
a.surv_rate.f <- function(temp.new) {
  surv_d <- c(2,52.33,76.77,66.33,5.87,10)
  temp_d <- c(5,18,23,28,33,35)
  model <- drm(surv_d ~ temp_d, fct = DRC.beta())
  a.surv.rate <- predict(model,data.frame(temp.v=temp.new))
  a.surv.rate <- 1-1/ifelse(a.surv.rate<1,1,a.surv.rate)
  return( a.surv.rate ) 
}
# plot(-15:45,a.surv_rate.f(-15:45),col="red",type="l")
# points(temp_d,1-1/surv_d,col="blue")
# cbind(-15:55,round(a.surv_rate.f(-15:55),2))

#emergence rate immature
#development: duration of development (number of days) from Tab. 2
i.emer_rate.f=function(temp.new){
  eme_d <- 1/c(35, 10.31, 4.19, 3.43,  3.00, 2.07, 1.82, 5, 25)
  temp_d <-  c(8,  13,    18,   23,    23,   28,   33,  35, 40)
  model <- drm(eme_d ~ temp_d, fct = DRC.beta())
  i.emer.pred <- predict(model,data.frame(temp.v=temp.new))
  return( i.emer.pred )
}
# plot(-15:50,i.emer_rate.f(-15:50),col="red",type="l")
# points(eme_d~temp_d,col="blue")
# cbind(-15:55,round(i.emer_rate.f(-15:55),2))

#immature survival rate. Marini et al. (2019): 80% at 8°C
i.surv_rate.f <- function(temp.new){
  surv_d <- c(0, 80, 97.8, 93.5, 95.6, 92.8, 84.7, 0)/100 
  temp_d <- c(4, 8, 13,   18,   23,   28,   33,   45)
  model <- drm(surv_d ~ temp_d, fct = DRC.beta())
  i.surv.pred <- predict(model,data.frame(temp.v=temp.new))
  return( i.surv.pred ) 
}
# plot(-15:50,i.surv_rate.f(-15:50),col="red",type="l")
# points(surv.v~temp.v,col="blue")
# cbind(-15:55,round(i.surv_rate.f(-15:55),2))

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
  hatc_d <- 1/c(100,2.45,1.35,1.07,1.08,2.04,100)/2
  temp_d <- c(7,8,13,23,28,33,35)
  model <- drm(hatc_d ~ temp_d, fct = DRC.beta())
  e.hatch.pred <- predict(model,data.frame(temp.v=temp.new))
  return( e.hatch.pred ) 
}
# plot(-10:50,e.hatch_rate.f(-10:50),col="red",type="l")
# points(hatc_d~temp_d,col="blue")
# cbind(-15:55,round(e.hatch_rate.f(-15:55),2))

#eggs survival rate.
e.surv_rate.f <- function(dt){
  ed_surv_bl=1
  a=0.98
  b=15.68
  c=-17.67
  d.surv.pred= ed_surv_bl*a*exp(-0.5*((dt-b)/c)^6)
  return( d.surv.pred )
}
# plot(-10:50,e.surv_rate.f(-10:50),col="red",type="l")
# points(hatc_d~temp_d,col="blue")
# cbind(-15:55,round(e.surv_rate.f(-15:55),2))

## Log-Normal probability density of short active dispersal (from DOI: 10.1002/ecs2.2977); from 0 to 600 m with resolution of 10 m.
f.adis.p <- dlnorm(seq(0,600,10), meanlog=4.95, sdlog=0.66)