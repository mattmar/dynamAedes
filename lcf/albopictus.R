#Temperature-dependent life cycle functions for `Aedes albopictus`

#emergence rate immature -> adult: from Tab. 1 in Delatte et al. (2009), column Pupae-adult
i.emer_rate.f <- function(dt){
  emerd <-c(0,0,83.3,89.9,93.8,90.0,8.3,0)
  st <- seq(5,40, by=5)
  model <- lm(emerd ~ poly(st,2))
  pred_emer_perc <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_perc[,1] <- ifelse(pred_emer_perc[,1]<0,0,pred_emer_perc[,1])
  pred_emer_rate <- pred_emer_perc[,1]/100
  return( pred_emer_rate )
}
# plot(5:40,i.emer_rate.f(5:40),col="red")
# points(emerd~st,col="blue")

## Immature daily survival rate at different temperature ##
## Data taken from:
i.mort_rate.f <- function(dt){
  a=0.977 
  b=21.8 
  c=-16.6
  pred_i_rate <- a*exp(-0.5*((dt-b)/c)^6)
  return( 1-(pred_i_rate) ) 
}
# plot(5:40,i.mort_rate.f(5:40),col="red")
# points(5:40,i.mort_rate.f(5:40),col="blue")

## Adult female daily mortality rate at different temperatures ##
## Data taken from:
a.mort_rate.f <- function(dt){
  a=0.677 
  b=20.9
  c=-13.2
  if( dt>0 ){
    pred_surv_rate <- a*exp(-0.5*((dt-b)/c)^6)*(dt^0.1)
  } else {
    pred_surv_rate <- a*exp(-0.5*((dt-b)/c)^6)
  }
  return( 1-(pred_surv_rate) ) 
}
# plot(5:40,a.mort_rate.f(5:40),col="red")
# points(5:40,a.mort_rate.f(5:40),col="blue")

## Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition. ##
# This process has a minimum duration of 5 days, which consiste of 24h to mature gonads, 3 days for gonadotrophic cycle.
## Data taken from: Delatte et al. JOURNAL OF MEDICAL ENTOMOLOGY (2009) 46:1; DOI: http://dx.doi.org/10.1603/033.046.0105 from Tab. 6 in Delatte et al. (2009)
a.gono_rate.f <- function(dt){
  gpd <-c(8.1,4.5,3.5,4.4)
  st <-c(20,25,30,35)
  model <- lm(gpd ~ poly(st,2))
  pred_gono_duration <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_gono_duration[,1] <- ifelse(pred_gono_duration[,1]<0,0,pred_gono_duration[,1])
  pred_gono_rate <- 1/pred_gono_duration[,1]
  return( pred_gono_rate ) 
}
# plot(0:40,a.gono_rate.f(0:40),col="red")
# points(0:40,a.gono_rate.f(0:40),col="blue")

## Oviposition rate, i.e., number of eggs laid per female/day at different temperature ##
## Data taken from: number of eggs from Tab. 6 in Delatte et al. (2009)
a.ovi_rate.f <- function(dt){
  ovir <-c(50.8, 65.3, 74.2, 48.7)
  st <-c(20,25,30,35)
  model <- lm(ovir ~ poly(st,2))
  pred_ovi_rate <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_ovi_rate[,1] <- ifelse(pred_ovi_rate[,1]<0,0,pred_ovi_rate[,1])
  return( pred_ovi_rate[,1] )
}
# plot(0:50,a.ovi_rate.f(0:50),col="red",ylim=c(0,100))
# points(0:50,a.ovi_rate.f(0:40),col="blue")

#Data taken from the Supplementary Materials available in : Metelmann et al. Journal of the Royal Society Interface (2019) 12:524; DOI: https://doi.org/10.6084/m9.figshare.c.4418279.v1
#eggs
e.surv_rate <- function(dt){
  a=0.955 
  b=18.8 
  c=-21.53
  pred_surv_rate= a*exp(-0.5*((dt-b)/c)^6)
  return(pred_surv_rate) 
}
# plot(0:50,e.surv_rate(0:50),col="red")
# points(0:50,e.surv_rate(0:40),col="blue")

# diapausing eggs (this is not yet exactely how they have it written)
ed.surv_rate <- function(dt){
  ed_surv_bl = 1
  a=0.93
  b=11.68
  c=-15.67
  pred_surv_rate= ed_surv_bl*a*exp(-0.5*((dt-b)/c)^6)
  return(pred_surv_rate)
}

# plot(0:50,e.surv_rate(0:50),col="red")
# points(0:50,ed.surv_rate(0:50),col="blue")
