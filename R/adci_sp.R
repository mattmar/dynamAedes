adci_sp=function (input_sim = NULL, coords = NULL, eval_date = NULL, stage = 1, breaks = c(0.25, 0.5, 0.75)){
  if (stage < 1 | stage > 3) {
    stop("stage can be 1 (egg), 2 (juvenile), 3 (adult) ...")
  }
  if (max(eval_date) > max(sapply(input_sim, length))) {
    stop("eval_date > than number of simulated days...")
  }
  if (all(unlist(lapply(input_sim, function(x) {
    sapply(x, length)
  })) == 4)) {
    stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes")
  }
  else {
    # subset life stage
    mylist = lapply(input_sim, function(x) {Map(function(z, y) z[y, ], x, stage)})
    mylist = lapply(mylist, function(x) { tmp = lapply(x, function(y) {  data.frame(y)})
            tmp = do.call(cbind, tmp)
            })
    # subset days
  myOut <- lapply(eval_date, function(z){

  mysublist=lapply(mylist, function(x){if(length(x)>= z){
    x1=x[z]
  }else{
    x1=cbind.data.frame(x, as.data.frame(matrix(0, nrow=nrow(x), ncol=z-ncol(x))))
  }
    return(x1)})
  
  # matrix operation
  my.mx <- lapply(mysublist, function(x) {
    if (any(class(x) == "data.frame")) 
      as.matrix(x)
    else x
  })
  my.out <- apply(simplify2array(my.mx), 1:2, quantile, 
                  prob = breaks)
  
  rate.sp <- data.frame(t(my.out[, , 1]))
  # rasterize
  rb <- lapply(1:length(breaks), function(y) {
      tmp = rate.sp[, y]
      tmp = as.data.frame(cbind(coords, tmp))
      names(tmp)[1:2] = c("X", "Y")
      tmp=rasterFromXYZ(tmp)
      names(tmp) = paste0("p", breaks[y])
      return(tmp)
        })
    outr=stack(rb)
})
  }
  names(myOut)=paste0("d_", eval_date)
  return(if (length(myOut) == 1) outr[[1]] else myOut)
}


