#' Function 1 to start Beta functions: Functions taken from https://github.com/OnofriAndreaPG/aomisc
#' @description Function 1 to start Beta functions
#' @keywords internal
#' @return No return value, called for side effects in .DRC.beta
.beta.fun <- function(X, b, d, Xb, Xo, Xc){
  .expr1 <-  (X - Xb)/(Xo - Xb)
  .expr2 <- (Xc - X)/(Xc - Xo)
  .expr3 <- (Xc - Xo)/(Xo - Xb)
  ifelse(X > Xb & X < Xc, d * (.expr1*.expr2^.expr3)^b, 0)
}

#' Function 2 to start Beta functions: Functions taken from https://github.com/OnofriAndreaPG/aomisc
#' @description Function 1 to start Beta functions
#' @keywords internal
#' @return No return value, called for side effects in drm
.DRC.beta <- function(){  
  fct <- function(x, parm) {
# function code here
.beta.fun(x, parm[,1], parm[,2], parm[,3], parm[,4], parm[,5])
}
ssfct <- function(data){
# Self-starting code here
x <- data[, 1]
y <- data[, 2]

d <- max(y)
Xo <- x[which.max(y)]
firstidx <- min( which(y !=0) )
Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
secidx <- max( which(y !=0) )
Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
c(1, d, Xb, Xo, Xc)

}
names <- c("b", "d", "Xb", "Xo", "Xc")
text <- "Beta function"

## Returning the function with self starter and names
returnList <- list(fct = fct, ssfct = ssfct, names = names, text = text)
class(returnList) <- "drcMean"
invisible(returnList)
}

#' Function 3 to start Beta functions: Functions taken from https://github.com/OnofriAndreaPG/aomisc
#' @description Function 1 to start Beta functions
#' @keywords internal
#' @return No return value, called for Beta function initiation in drm
.beta.init <- function(mCall, LHS, data) {
  xy <- sortedXyData(mCall[["X"]], LHS, data)
  x <-  xy[, "x"]; y <- xy[, "y"]

#Self starting code ##############
d <- max(y)
Xo <- x[which.max(y)]
firstidx <- min( which(y !=0) )
Xb <- ifelse(firstidx == 1,  x[1], (x[firstidx] + x[(firstidx - 1)])/2)
secidx <- max( which(y !=0) )
Xc <- ifelse(secidx == length(y),  x[length(x)], (x[secidx] + x[(secidx + 1)])/2)
start <- c(1, d, Xb, Xo, Xc)
names(start) <- mCall[c("b", "d", "Xb", "Xo", "Xc")]
start
}
#NLS.beta <- selfStart(.beta.fun, .beta.init, parameters=c("b", "d", "Xb", "Xo", "Xc"))

#' Function 4 to avoid NAs in the predicted rates
#' @description function to clean predicted rates
#' @keywords internal
#' @return vector of rates.
#' @noRd
sanitize_output <- function(x, is_prob = TRUE, min_val = 1e-6, max_val = 1) {
  # Replace NA and NaN with min_val (safe default)
  x[is.na(x)] <- min_val
  
  # Replace -Inf with min_val (very unsafe = no survival)
  x[is.infinite(x) & x < 0] <- min_val
  
  # Replace +Inf with max_val (very safe = full survival)
  x[is.infinite(x) & x > 0] <- max_val
  
  # Clamp only if outside expected bounds
  if (is_prob) {
    x <- pmin(pmax(x, 0), 1)
  } else {
    x <- pmax(x, 0)
  }
  
  return(x)
}

#' Daily rate of gonotrophic cycle
#' @description Daily rate of gonotrophic cycle, i.e., daily rate of adult females completing the cycle from bloodmeal to the first day of oviposition
#' @keywords internal
#' @return vector of rates.
.a.gono_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    gono_d <- 1/c(30,15,7,5,4,4,5,5,6,7,35,100)
    temp_d <- c(15,17,20,22,26,28,31,32,33,33,40, 45)
    model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
    a.gono.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
    }else if(sp=="albopictus") {
      gono_d <- 1/c(8.1, 4.5, 3.5, 4, 4.4, 100)
      temp_d <-   c(20,  25,  30,  32.5, 35, 45)
      model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
      a.gono.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
      }else if(sp=="koreicus") {
        gono_d <- 1/c(14.75, 11.5, 9.21, 10, 10.81, 100)
        temp_d <- c(  18,    23,   23,    25, 28,    33)
        model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
        a.gono.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
        }  else if(sp=="japonicus") {
          gono_d <- 1/c(14.75, 11.5, 9.21, 10, 10.81, 100)
          temp_d <- c(  18,    23,   23,    25, 28,    33)
          model <- drm(gono_d ~ temp_d, fct = .DRC.beta())
          a.gono.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
          }  else(stop("Species not supported."))
          return(a.gono.pred)
        }

#' Oviposition rate
#' @description Oviposition rate, i.e., number of eggs laid per female/day at different temperature
#' @keywords internal
#' @return vector of rates.
.a.ovi_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    ovi_n <- c( 0, 2.5, 3.37, 4.6,  6.98, 7.6,  9.58, 7.28,11.22,7.27,0) * 10
    temp_n <- c(0, 10,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,34)
    model <- lm(ovi_n~poly(temp_n,3))
    a.ovi.pred <- predict(model, newdata = data.frame(temp_n=temp.new), response=TRUE)
    a.ovi.pred <- ifelse(a.ovi.pred<0, 0, a.ovi.pred)
    a.ovi.pred[which(temp.new<4|temp.new>45)] <- 0
    }else if(sp=="albopictus") {
      ovi_n <-  c(0, 50.8, 65.3, 74.2, 48.7, 0)
      temp_n <- c(5, 20,   25,   30,   35, 45)
      model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
      a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
      }else if(sp=="koreicus") {
        ovi_n <- c(0, 5, 20, 25, 30, 38, 40, 40, 38, 20, 10, 10, 0)
        temp_n <-c(8, 10,12, 15,17, 20, 23,  25, 27, 30, 33, 35,37)
        model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
        a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
        }else if(sp=="japonicus") {
          ovi_n <- c(0, 108.2, 111.6, 106.8, 112.2, 97.1, 99.1, 94.5, 80.6, 82.1, 71.6, 67.4, 68.4, 55.0, 47.4,0)
          temp_n <-c(5, 10, 12, 14, 15, 17, 19, 20, 23, 25, 26, 27, 28, 29, 31, 40)
          model <- drm(ovi_n ~ temp_n, fct = .DRC.beta())
          a.ovi.pred <- predict(model,data.frame(temp.v=temp.new))
          }else(stop("Species not supported."))
          return( a.ovi.pred )
        }

#' Adult female daily mortality rate
#' @description Adult female daily mortality rate at different temperatures
#' @keywords internal
#' @return vector of rates.
.a.surv_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    surv_d <- c(3,13.18,10.91,27.71,30.62,23.72,26.90,32.87,36.91,22.77,29.26,22.53,10.07,5)
    temp_d <- c(5,10.54,10.76,15.30,16.52,20.05,21.79,25.64,27.64,31.33,31.65,32.55,33.41,36)
    model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
    a.surv.rate <- predict(model,data.frame(temp.v=temp.new))
    a.surv.rate <- 1-1/ifelse(a.surv.rate<1,1,a.surv.rate)
    a.surv.rate<-sanitize_output(a.surv.rate)
    }else if(sp=="albopictus") {
      a=0.677; b=20.9; c=-13.2
      a.surv.rate <- sapply(temp.new, function(x) {
        if( x>0 ){
          a*exp(-0.5*((x-b)/c)^6)*(x^0.1)
          } else {
            a*exp(-0.5*((x-b)/c)^6)
          }
          })
      }else if(sp=="koreicus" | sp=="japonicus") {
        surv_d <- 1-(1/c(1,2,52.33,46.77,66.33,5.87, 3, 1))
        temp_d <- c(0,5, 18,23,28,33, 35, 38)
        model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
        a.surv.rate <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
#koreicus adults longevity (Tab. 3, Marini et al 2019)
temp_dev <- c(1,5,18,23,28,33) #I've added arbitrary adult longevity at 5°C equal to 30 days, and adult longevity at 1°C equal to 1 days 
obs_dev <- c(1, 30, 52.33, 46.77, 66.33, 5.87)
m <- lm(obs_dev~poly(temp_dev,2))
dev_time <- predict(m, newdata = data.frame(temp_dev=temp.new), response=TRUE)
dev_time <- ifelse(dev_time<=1, 1, dev_time)
#correction with development time to get daily survival rate 
a.surv.rate <- round(a.surv.rate^(1/dev_time),3)
a.surv.rate <- sanitize_output(a.surv.rate)
}else if(sp=="japonicus") {
  surv_d <- c(0, 99.0, 97.2, 91.8, 98.8, 91.0, 81.6, 93.6, 68.9, 71.4, 71.4, 60.7, 68.9, 34.0, 40.8,0 )/100
  temp_n <- c(-5, 10, 12, 14, 15, 17, 19, 20, 23, 25, 26, 27, 28, 29, 31, 40)
  model <- drm(surv_d ~ temp_n, fct = .DRC.beta())
  a.surv.rate <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
}
else(stop("Species not supported."))
return( a.surv.rate )
}

#' Log-Normal probability density of short active dispersal
#' @description Log-Normal probability density of short active dispersal from Marcantonio et al. 2019; Marini et al. 2019. 0 to 600 m with resolution of 10 m. 
#' @keywords internal
#' @return vector of rates.
.a.a_disp.f <- function(sp, max.a.disp, disp.bins) {
  if(sp=="aegypti"){
    dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.95, sdlog=0.66)
    }else if(sp=="albopictus") {
      dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.54, sdlog=0.58)
      }else if(sp=="koreicus") {
dispk <- dlnorm(seq(0,max.a.disp,disp.bins), meanlog=4.54, sdlog=0.58) #check paper svizzera
}else(stop("Species not supported."))
return( dispk )
}

#' Immature daily emergence rate
#' @description Immature daily emergence rate at different temperature
#' @keywords internal
#' @return vector of rates.
.i.emer_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    eme_d <- 1/(c(500,61.7,39.7,84.4,10.0,9.2,8.4,6.3,7.4,5.1,8.1,6.3,50)*0.2980392)
    temp_d <- c(12,14.74,14.84,14.92,26.56,26.84,26.85,30.83,31.61,34.95,36.47,36.55,42)
    model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
    i.emer.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
    }else if(sp== "albopictus") {
      eme_d <- 1/c(100,8.7,4.1,2.7,1.9,1.7,5)
      temp_d <- c(5,15,20,25,30,35,38)
      model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
      i.emer.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
      }else if(sp=="koreicus") {
        eme_d <- 1/c(35, 10.31, 4.19, 3.43,  3.00, 2.07, 1.82, 5, 25)
        temp_d <-  c(8,  13,    18,   23,    23,   28,   33,  35, 40)
        model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
        i.emer.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
        }else if(sp=="japonicus") {
# This refers to the total % of emergency, it must be corrected with the durations to get a daily rate
eme_d= c( 0.07, 0.1, 0.15, 0.17, 0.21, 0.32, 0.29, 0.35, 0.43, 0.43, 0.44, 0.45, 0.5, 0.67, 0.4,0.2, 0)
temp_d <- c(10, 12, 14, 15, 17, 19, 20, 23, 25, 26, 27, 28, 29, 31, 33, 35,40)
model <- drm(eme_d ~ temp_d, fct = .DRC.beta())
i.emer.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
}else(stop("Species not supported."))
return( i.emer.pred )
}

#' Immature daily survival rat
#' @description Immature daily survival rate at different temperature
#' @keywords internal
#' @return vector of rates.
.i.surv_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    surv_d <- 1-(1/ c(1,6.2,10.3,12,2,4.5,2.1,5.7,4.8,47.9,42.7,55.3,27.2,48.5,30.2,14.7,15.8,15.7,9.5,16.8,8.6,14.2,1))
    temp_d <- c(2,10,10,10,10,10.38,10.45,10.45,10,14.74,14.84,14.92,18.86,19.04,19.18,26.56,26.84,26.85,30.83,31.61,34.95,36.47,39)
    model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
    # pred <- suppressWarnings(predict(model, newdata = data.frame(temp.v = temp.new)))
        # Replace any NA/NaN/negative/zero predictions with small fallback value
    # pred[is.na(pred) | pred <= 0 | !is.finite(pred)] <- 1e-6
    i.surv.pred <- sanitize_output(predict(model, newdata = data.frame(temp.v = temp.new)))
    }else if(sp=="albopictus") {
      a=0.977
      b=20.8 
      c=-12.6
      i.surv.pred <- a*exp(-0.5*((temp.new-b)/c)^6)
      }else if(sp=="koreicus") {
#surv_d <- c(0, 50, 80, 97.8, 93.5, 95.6, 92.8, 84.7, 0)/100 # patch
#temp_d <- c(-5, 0, 8, 13,   18,   23,   28,   33,   37) #patch 
surv_d <- c(0, 80, 97.8, 93.5, 95.6, 92.8, 84.7, 0)/100 
temp_d <- c(-5, 8, 13,   18,   23,   28,   33,   37)
model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
i.surv.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
#koreicus adults longevity (Tab. 2, Marini et al 2019)
temp_dev <- c(13, 18,23,28,33)
obs_dev <- c(9.47,5.71,3.66,2.72,2.78)
m <- lm(obs_dev~poly(temp_dev,2))
dev_time <- predict(m, newdata = data.frame(temp_dev=temp.new), response=TRUE)
i.surv.pred <- sanitize_output(round(i.surv.pred^(1/dev_time),3))
} else if(sp=="japonicus") {
  surv_d <- c(0, 0, 92.5, 78.77, 90.59, 92.19, 90.32, 83.86, 94.08, 91.35, 75.5, 90.63, 96.95, 90.83, 53, 38.65, 0)/100
  temp_d <- c(0, 5, 10, 12, 14, 15, 17, 19, 20, 21, 23, 25, 26, 27, 29, 31,45)
  model <- drm(surv_d ~ temp_d, fct = .DRC.beta())
  i.surv.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
  }else(stop("Species not supported."))
  return( i.surv.pred )
}

#' Immature Density-dependent mortality rate
#' @description Immature Density-dependent mortality
#' @keywords internal
#' @return vector of rates.
.i.ddmort_rate.f <- function(dens.new) {
  i.dens.v = c(87.72,561.40,1175.44,1280.70,1491.23,1675.44,1982.46,2350.88,2850.88,3122.81,3236.84,3307.02,3359.65,3456.14,3570.18,3640.35,3666.67,3771.93,3877.19,3982.46)
  i.dsurv_prop.v = c(0.95,0.43,0.57,0.42,0.49,0.35,0.25,0.17,0.13,0.05,0.2,0.27,0.11, 0.11,0.06,0.04,0.09,0.13,0.07, 0.14)
  i.dsur_dur.v = c(7.46,9.96,28.09,17.46,29.12,38.31,42.22,40.95,44.31,42.34,20.91,17.43,
    37.12,29.73,40,43,22,51,50,31)
# Transform survival proportion to multidays mortality rate
i.dmort_rate.v = -log(i.dsurv_prop.v)
# Transform multiday mortality rate to daily mortality rate
i.dmort_rate.v = i.dmort_rate.v/i.dsur_dur.v
# Fit a model for daily mortality rate using just density as a predictor
i.dmort.lm <- lm(log(i.dmort_rate.v) ~ i.dens.v)
i.dmort.pred <- predict(i.dmort.lm, dens.new)
return( i.dmort.pred )
}

#' Egg hatching rate
#' @description Egg hatching rate: this rate decides embryonated eggs which hatch or stay
#' @keywords internal
#' @return vector of rates.
.e.hatch_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    hatc_d <- c(0,0, 0.025,0.98, 0.99, 0.73, 0.30, 0.016, 0)
    temp_d <- c(7,12,19,   23,   24.5, 26.5, 29.5, 40,    45)
    model <- lm(hatc_d ~ poly(temp_d,5))
    e.hatch.pred <- sanitize_output(predict(model,data.frame(temp_d=temp.new), response=TRUE))
    e.hatch.pred[which(temp.new<=15 | temp.new>=39)] <- 0
    }else if(sp=="albopictus") {
      hatc_d <-c(2,4.4,4.0,50,100,60,51.4,10.0,0)/100 
      temp_d <-seq(0,40,by=5)
      model <- lm(hatc_d ~ poly(temp_d,4))
      e.hatch.pred <- sanitize_output(predict(model,data.frame(temp_d=temp.new), response=TRUE)) 
      e.hatch.pred <- ifelse(e.hatch.pred<0, 0, e.hatch.pred)
      e.hatch.pred[which(temp.new<4|temp.new>40)] <- 0
      }else if(sp=="koreicus") {
        hatc_d <- (c(0, 7.25, 60.50, 43.75, 31.00, 27.25,10,0)/100)/0.636
        temp_d <- c(0,8,13,23,28,33,36, 36.5)
        model <- lm(hatc_d ~ poly(temp_d,3))
        e.hatch.pred <- sanitize_output(predict(model,data.frame(temp_d=temp.new), response=TRUE))
        e.hatch.pred <- ifelse(e.hatch.pred<0, 0, e.hatch.pred)
        e.hatch.pred[which(temp.new<4)] <- 0
        }else if(sp=="japonicus") {
          hatc_d <- c(0.4075, 0.8275, 0.9175, 0.89,0)  
          temp_d <- c(0,10,20,30,35)
          model <- drm(hatc_d ~ temp_d, fct = .DRC.beta())
          e.hatch.pred <- sanitize_output(predict(model,data.frame(temp.v=temp.new)))
          }else(stop("Species not supported."))
          return( e.hatch.pred )
        }

#' Egg daily survival rate
#' @description Egg daily survival rate at different temperature
#' @keywords internal
#' @return vector of rates.
.e.surv_rate.f <- function(temp.new, sp) {
  if(sp=="aegypti") {
    surv_r <- c(0,0,0,0,0.2,0.3,0.5,0.7,0.78,0.81,0.88,0.95,1.00,0.93,0.93,0.83,0.90,0.48,0)
    temp_r <- c(-17,-15,-12,-10,-7,-5,-2,0,15.6,16,21,22,25.0,26.7,28.0,31.0,32.0,35.0,41.0)
    model <- lm(surv_r ~ poly(temp_r,5))
    e.surv.pred <- sanitize_output(predict(model,data.frame(temp_r=temp.new)))
    e.surv.pred[temp.new >= 41 | temp.new <= -10] <- 0
    }else if(sp=="albopictus") {
      a= 0.955
      b= 16
      c= -17 
      e.surv.pred=a*exp(-0.5*((temp.new-b)/c)^6)
      }else if(sp=="koreicus") {
        ed_surv_bl=1
        a= 0.98
        b= 15.8
        c= -15.8
        e.surv.pred=ed_surv_bl*a*exp(-0.5*((temp.new-b)/c)^6)
        }else if(sp=="japonicus") {
          surv_r <-  c(0, 0.44, 0.79, 0.9017, 0.8817, 0)
          temp_r <-c(-15,0,10,20,30,35)
          model <- drm(surv_r ~ temp_r, fct = .DRC.beta())
          e.surv.pred <- predict(model,data.frame(temp.v=temp.new))
          }else(stop("Species not supported."))
          return( e.surv.pred )
        }

#' Diapause eggs daily survival rate
#' @description Diapause eggs daily survival rate at different temperature
#' @keywords internal
#' @return vector of rates.
.d.surv_rate.f <- function(temp.new, sp){
  if(sp=="albopictus") {
    ed_surv_bl=1
    a=0.955
    b=11.68
    c=-15.67
    d.surv.pred=ed_surv_bl*a*exp(-0.5*((temp.new-b)/c)^6)
    }else if(sp=="koreicus" | sp=="japonicus") {
      d.surv.pred=rep(0.999,length(temp.new))
      }else(stop("Species not supported."))
      return( d.surv.pred )
    }

#' Allocation of diapause eggs
#' @description Allocation of diapause eggs dependent on photoperiod
#' @keywords internal
#' @return vector of rates.
.e.dia_rate.f <- function(photo.new, sp) {
  if(sp=="albopictus") {
    e.diap.pred <-1/(1+exp(3.04*(photo.new-12.62)))
    }else if(sp=="koreicus" | sp=="japonicus") {
      e.diap.pred <-1/(1+exp(3.04*(photo.new-12.97)))
      }  else(stop("Species not supported."))
      return(e.diap.pred)
    }

#' Distance moved by mosquito populations
#' @description Return quantiles of distance moved by mosquito populations
#' @keywords internal
#' @return vector of quantiles of distances.
.returndis <- function(distl=NA, days=NA, breaks=breaks) {
  out <- do.call(rbind.data.frame, 
    lapply(days, function(x) {
      round(quantile(unlist(sapply(1:length(distl),
        function(y) {mean(if(x<=length(distl[[y]])) as.integer(distl[[y]][[x]]) else NA, na.rm=T)})), probs=breaks, na.rm=T),1)
      }))
  names(out) <- breaks
  return(out)
}

#' Euclidean distance
#' @description Return Euclidean distance
#' @keywords internal
#' @return float. An euclidean distance.
.euc <- function(xs, ys) { sqrt((xs[1]-xs[2])^2 + (ys[1]-ys[2])^2) }

#' Euclidean distance from a pair of coordinates
#' @description Returns Euclidean distance from a pair of coordinates
#' @keywords internal
#' @return float. Vector of Euclidean distances.
.meuc <- function(c1, c2, coords) {
  c3 <- coords[unlist(c1[,1])[1],]
  outd <- list(NA)
  for ( rw in 1:length(c2) ) {
    outd[[rw]] <- sapply(unlist(c2[rw]), function(x) {
      .euc(c(as.numeric(c3[1]), as.numeric(coords[x,1])), c(as.numeric(c3[2]), as.numeric(coords[x,2])))
      })
  }
  return(outd)
}

#' Safer version of sample
#' @description sample
#' @keywords internal
#' @return  a vector of length n with elements drawn from either x or from the integers 1:x
.resample <- function(x, ...) x[sample.int(length(x), ...)]

#' Date format check
#' @description Check the proper format of a date object
#' @keywords internal
#' @return Logical
#' 
# Check dayspan as well as date format function
check_date_format <- function(date) {
  if (!grepl("^\\d{4}-\\d{2}-\\d{2}$", date)) {
    stop("Dates in the wrong format: change them to '%Y-%m-%d'.")
  }
}