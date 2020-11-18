setwd('/home/matteo/own_data/PoD/topics/aedes_genmod/lcf')
library(epiR)
temp.v <- seq(-15,45,0.1)
dev.off()

png("~/Aedes_lifecycle_functions_beta.png",heigh=480*2,width=480*5.5,res=100)
par(pch=16,cex=0.5,cex.lab=2,cex.axis=1.5,cex.main=2,mfrow=c(2,4))
# Gonotrophic rate
source("./aegypti_beta.r")
plot(temp.v,a.gono_rate.f(temp.v),col="red",type="l",ylim=c(0,0.3),ylab="Gonotrophic rate",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,a.gono_rate.f(temp.v),col="blue")
source("./koreicus_beta.r")
lines(temp.v,a.gono_rate.f(temp.v),col="green")
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily gonotrophic rate")

# Daily fertility rate
source("./aegypti_beta.r")
plot(temp.v,a.ovi_rate.f(temp.v),col="red",type="l",ylim=c(0,110),ylab="Number of eggs/female",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,a.ovi_rate.f(temp.v),col="blue",ylim=c(0,1))
source("./koreicus_beta.r")
lines(temp.v,a.ovi_rate.f(temp.v),col="green",ylim=c(0,1))
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily fertility rate")

# Daily survival rate
source("./aegypti_beta.r")
plot(temp.v,a.surv_rate.f(temp.v),col="red",type="l",ylim=c(0,1),ylab="Adult survival",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,a.surv_rate.f(temp.v),col="blue",ylim=c(0,1))
source("./koreicus_beta.r")
lines(temp.v,a.surv_rate.f(temp.v),col="green",ylim=c(0,1))
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily adult survival rate")

# Daily immature emergency rate
source("./aegypti_beta.r")
plot(temp.v,i.emer_rate.f(temp.v),col="red",type="l",ylim=c(0,0.7),ylab="Immature emergence",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,i.emer_rate.f(temp.v),col="blue")
source("./koreicus_beta.r")
lines(temp.v,i.emer_rate.f(temp.v),col="green")
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily immature emergency rate")

# Daily immature density-independent survival rate
source("./aegypti_beta.r")
plot(temp.v,i.surv_rate.f(temp.v),col="red",type="l",ylim=c(0,1),ylab="Immature survival",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,i.surv_rate.f(temp.v),col="blue")
source("./koreicus_beta.r")
lines(temp.v,i.surv_rate.f(temp.v),col="green")
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily immature survival rate")

# Daily immature density-dependant survival rate
source("./aegypti.r")
plot(seq(1,50000,100)/4, 1-(1-exp(-exp(predict(i.ddmort_rate.m,data.frame(i.dens.v=seq(1,50000,100)))))) ,type="l", ylab="Survival daily rate",xlab="Abundance per L of water")
legend("topright",y=0.5,c("aegypti","albopictus","koreicus"),cex=2,col=c("black"),lt=1,bty="n")
title("Daily immature\ndensity-dependent survival rate")

# Daily hatching rate
source("./aegypti_beta.r")
plot(temp.v,e.hatch_rate.f(temp.v),col="red",type="l",ylim=c(0,1),ylab="Hatching rate",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,e.hatch_rate.f(temp.v),col="blue")
source("./koreicus_beta.r")
lines(temp.v,e.hatch_rate.f(temp.v),col="green")
legend("topleft",c("aegypti","albopictus","koreicus"),cex=2,col=c("red","blue","green"),lt=1,bty="n")
title("Daily egg hatching rate")

# Daily egg survival rate
source("./aegypti_beta.r")
plot(temp.v,e.surv_rate.f(temp.v),col="red",type="l",ylim=c(0,1),ylab="Survival rate",xlab="Temperature °C")
source("./albopictus_beta.r")
lines(temp.v,e.surv_rate.f(temp.v),col="blue")
lines(temp.v,d.surv_rate.f(temp.v),col="cyan3")
source("./koreicus_beta.r")
lines(temp.v,e.surv_rate.f(temp.v),col="green")
legend(x=6,y=0.32,c("aegypti","albopictus","albopictus dia","koreicus"),cex=2,col=c("red","blue","cyan3","green"),lty=1,bty="n")
title("Daily egg survival rate")

dev.off()

