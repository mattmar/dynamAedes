#' Temporal downscaling of entomological surveillance observations
#'
#' @param mydf A data.frame.
#' @param date.field character, name of field containing dates.
#' @param value.field integer, name of field containing the number of individuals sampled.
#' @param counter.field integer, name of the field containing the number of days/weeks between each sampling. 
#' This is optional, if NULL then the function will compute the value assuming the trap was active during the whole 
#' period between two consecutive samplings.
#' @param seed integer, define the seed for the binomial draws, default \code{seed = 123}.
#' @return Returns a data.frame with the "adjusted value", i.e. the observation spread over the period of activity of the trap
#' @author Daniele Da Re \email{dare.daniele@gmail.com}, Giovanni Marini \email{dare.daniele@gmail.com}
#' @export

spreader<-function(mydf = NULL, date.field=NULL, value.field=NULL, counter.field=NULL, seed=123){
  set.seed(seed)
  mydf$month <- as.numeric(format((as.Date(mydf[[date.field]])), "%m"))
  mydf$value_adj <- NA 
  
  #winter correction
  suppressWarnings({
  mydf$value_adj[mydf$month==1] <- sum(mydf[which(mydf$month==1),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==1),][["value"]])
  mydf[which(mydf$month==1 & is.na(mydf$value_adj)),][["value_adj"]] <- 0 
  mydf$value_adj[mydf$month==2] <- sum(mydf[which(mydf$month==2),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==2),][["value"]])
  mydf[which(mydf$month==2 & is.na(mydf$value_adj)), ][["value_adj"]] <- 0 
  })
  
  #if statement to end function in case all the values are 0 or NAs
  if(sum(mydf[, value.field], na.rm=TRUE) <=10) {
    mydf$counter <- NA
    # complete dataset with observations at the end of the year
    suppressWarnings({
      mydf[which(mydf$month==10 & is.na(mydf$value_adj)),][["value_adj"]] <- 0 #to avoid NA in Oct due sampling activities between Oct and Nov
      mydf[which(mydf$month==11),][["value_adj"]] <- sum(mydf[which(mydf$month==11),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==11),][["value"]])
      mydf[which(mydf$month==12),][["value_adj"]] <- sum(mydf[which(mydf$month==12),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==12),][["value"]])
      mydf[which(mydf$month>=11),][["value_adj"]] <- ifelse(is.na(mydf[which(mydf$month>=11),][["value"]]) & is.na(mydf[which(mydf$month>=11),][["value_adj"]]) , 0, mydf[which(mydf$month>=11),][[ "value_adj"]])
    })  
  } else {
    # section of the function creating the counter (dataset without known revisiting period) 
    if(is.null(counter.field)){
      mydf$counter <- NA
      mydf$counter[min(which(!is.na(mydf[,value.field]) & mydf$month>=3 ))-1] <- 1 
      for(i in min(which(!is.na(mydf[,value.field]) & mydf$month>=3 )):nrow(mydf)){ 
        if (mydf$month[i] < 11){
          if (is.na(mydf[i,][[value.field]])){
            mydf$counter[i] <- mydf$counter[i-1]+1
          } else {
            mydf$counter[i] <- mydf$counter[i-1]+1
            activ.period <- rev(seq(1:(mydf$counter[i])))
            egg.size <- mydf[i,][[value.field]]
            sample.prob <- length(activ.period)
            
            for(j in activ.period){
              mydf[["value_adj"]][i-j+1] <- rbinom(size=egg.size, n=1, prob=1/sample.prob)   
              egg.size <- egg.size-mydf[["value_adj"]][i-j+1]
              sample.prob <- sample.prob-1
            }
          mydf$counter[i] <- 0
          }
        }
      }
    # section of the function working on the counter provided by the user
    } else {
      mydf$counter <- mydf[[counter.field]]
      for(i in min(which(!is.na(mydf[,value.field]) & mydf$month>=3 )):nrow(mydf)){ 
          if (mydf$month[i]<11 & !is.na(mydf[[value.field]][i])){
          activ.period <- rev(seq(1:(mydf$counter[i])))
          egg.size <- mydf[[value.field]][i]
          sample.prob <- length(activ.period)
          
          for(j in activ.period){
            mydf[["value_adj"]][i-j+1] <- rbinom(size=egg.size, n=1, prob=1/sample.prob)   
            egg.size <- egg.size-mydf[["value_adj"]][i-j+1]
            sample.prob <- sample.prob-1
          }
        } 
      }
    }
    # complete dataset with observations at the end of the year
    suppressWarnings({
      mydf[which(mydf$month==10 & is.na(mydf$value_adj)),][["value_adj"]] <- 0 #to avoid NA in Oct due sampling activities between Oct and Nov
      mydf[which(mydf$month==11),][["value_adj"]] <- sum(mydf[which(mydf$month==11),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==11),][["value"]])
      mydf[which(mydf$month==12),][["value_adj"]] <- sum(mydf[which(mydf$month==12),][["value"]], na.rm = TRUE)/length(mydf[which(mydf$month==12),][["value"]])
      mydf[which(mydf$month>=11),][["value_adj"]] <- ifelse(is.na(mydf[which(mydf$month>=11),][["value"]]) & is.na(mydf[which(mydf$month>=11),][["value_adj"]]) , 0, mydf[which(mydf$month>=11),][[ "value_adj"]])
    })  
    }
  mydf<- mydf[,!(names(mydf) %in% c("counter", "month"))  ]
  return(mydf)
}