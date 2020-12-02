#Aedes koreicus survival rates
#Daniele Da Re, Matteo Marcantonio 
#Data taken from: Marini et al. Parasites Vectors (2019) 12:524; DOI: https://doi.org/10.1186/s13071-019-3772-5
dt=seq(-10,35, 0.5)

#adults surv rate. We used the function described as F3 and parameters value present in Tab. 5 
a.koreicus_surv_rate=function(dt){
  a=0.01
  b=5.6 *10^-9
  c=0.52
 pred_surv_rate= a+b*exp(c*dt)
 pred_surv_rate=1-pred_surv_rate
 return(pred_surv_rate) 
}
plot(dt,a.koreicus_surv_rate(dt) )

#check for temperatures below 15 degress: Kim et al + expert based

#immature survival rate. We considered the Mean values of Instars and Pupae from Tab.1
i.koreicus_surv_rate=function(dt){
  emerd <-c(0.00 , 0.00 ,   0 ,  NA , 97.772,93.462,95.584,92.778,84.658) 
  st <-c(-10,   2,   4,   8,  13,  18,  23,  28,  33)
  model <- lm(emerd ~ poly(st,3))
  pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])
  pred_emer_rate <-pred_emer_time[,1]/100
  # plot(dt, pred_emer_rate)
  return( pred_emer_rate ) 
}
plot(dt,i.koreicus_surv_rate(dt) )

#eggs survival rate. We used the function described as F3 and parameters value present in Tab. 5 
e.koreicus_surv_rate=function(dt){
  a=0.77 
  b=7308.51 
  c=-0.93
  pred_surv_rate= a+b*exp(c*dt)
  pred_surv_rate= 1-pred_surv_rate
  return(pred_surv_rate) 
}
plot(dt,e.koreicus_surv_rate(dt) )

#eggs hatching rate
e.koreicus_hatc_rate=function(dt){
  emerd <-c(0.00 , 0.00 ,   NA , 7.25 ,50.50, 52.13, 53.75, 51.00, 57.25)
  st <-c(-10,   2,   4,   8,  13,  18,  23,  28,  33)
  model <- lm(emerd ~ poly(st,4))
  pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])
  pred_emer_rate <-pred_emer_time[,1]/100
  # plot(dt, pred_emer_rate)
  # points(st, emerd/100, bg="red", pch=21)
  return( pred_emer_rate ) 
}

plot(dt,e.koreicus_hatc_rate(dt) )
points(st, emerd/100, bg="red", pch=21)


#gonotrophic cycle (Tab. 3). Marini et al. (2019) assumed a gonotrophic cycle to last 11.6 days, the average computed with the experiments between 18 °C and 28 °C
koreicus_gonotrophic=function(dt){
  emerd <-c(14.75 , 9.21, 10.81 , 0)
  st <-c(18,23,28,33)
  # plot(st, emerd)
  model <- lm(emerd ~ poly(st,2))
  pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])
  pred_emer_rate <-pred_emer_time[,1]
  # plot(dt, pred_emer_rate)
  return( pred_emer_rate ) 
}
plot(dt,koreicus_gonotrophic(dt), ylim = c(0,10) )

#oviposition: In Marini et al. (2019), the average number of eggs laid in one ovipositionwas set to 100 accordingly to:# 
# Ciocchetta S, Darbro JM, Frentiu FD, Montarsi F, Capelli G, Aaskov JG,  et al. Laboratory colonization of the European invasive mosquito Aedes
# (Finlaya) koreicus. Parasit Vectors. 2017;10:74.



#development: duration of development (number of days) from Tab. 2

e.koreicus_dev=function(dt){
  emerd <-c(2.45,1.35, NA, 1.07,1.08,1.04)
  st <-c(8,13,18,23,28,33)
  # plot(st, emerd)
  model <- lm(emerd ~ poly(st,3))
  pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])
  pred_emer_rate <-pred_emer_time[,1]
  # plot(dt, pred_emer_rate)
  # points(st, emerd, bg="red", pch=21)
  return( pred_emer_rate )
}
plot(dt,e.koreicus_dev(dt) )

#development: duration of development (number of days) from Tab. 2
#Immature: average for colums L1 to Pupae
i.koreicus_dev=function(dt){
  emerd <-c(NA, 9.062, 4.542, 2.95, 2.42,2.486 )
  st <-c(8,13,18,23,28,33)
  plot(st, emerd)
  model <- lm(emerd ~ poly(st,3))
  pred_emer_time <- predict(model,newdata=data.frame(st=dt),interval='confidence',level=0.95)
  pred_emer_time[,1] <- ifelse(pred_emer_time[,1]<0,0,pred_emer_time[,1])
  pred_emer_rate <-pred_emer_time[,1]
  plot(dt, pred_emer_rate)
  points(st, emerd, bg="red", pch=21)
  return( pred_emer_rate )
}
plot(dt,i.koreicus_dev(dt) )


#hatcing rate, Tab 4
# e.koreicus_hatch_rate2=function(dt){
#   a=0.15 
#   b=1.28 
#   pred_surv_rate= exp(a-(b*dt))/(1+ exp(a-(b*dt)))
#   pred_surv_rate= 1-pred_surv_rate
#   return(pred_surv_rate) 
# }
# plot(dt,e.koreicus_hatch_rate2(dt) )


