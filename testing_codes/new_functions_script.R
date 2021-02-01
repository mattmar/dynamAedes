## Derive prop of a successfull introduction at the end of the simulated period
prob_succ_intro=function(input_sim=NULL, eval.date=NULL, st= NULL){ 
  #st can be: 1 = eggs, 2= immatures; 3=adults; 4=diapaused eggs
  
  pe_out=list()
  for(i in 1:length(eval.date)){
     pe=sum(unlist(lapply(input_sim, function(x) {
      pe <- length(which(sum(x[eval.date[i]][[st]])>0)); # when st >0 at a given date
    }))) / length(input_sim)
    pe_out[[i]]=pe
  }
  
  pe_out=data.frame(do.call(rbind,pe_out))
  names(pe_out)="succ_intro"
  rowid=paste0("day_after_intro_", eval.date)
  pe_out=data.frame(cbind("dai"=rowid, "succ_intro"=pe_out))
  pe_out$stage=rep(st, nrow(pe_out))
  return(pe_out)
  
}

## SPATIAL: Derive prop of a successfull introduction at the end of the simulated period 
prob_succ_intro_spatial = function(coords=NULL, input_sim=NULL ,eval.date=NULL){
    mylist=lapply(input_sim, function(x){
    mydf=lapply(x, function(y) {data.frame(y); data.frame("tot_individuals"= colSums(y))});
    mydf=do.call(cbind, mydf);
  })
  
  #add days and convert to 1 all pixel having a number of individuals >1
  col_names=paste0("d_", 1:ncol(mylist[[1]]))
  mylist=lapply(mylist, setNames, col_names)
  mylist=lapply(mylist, function(x) {x[x>0] <- 1; return(x)})
  
  #subset
  my_names=col_names[eval.date]
  mylist = lapply(mylist, "[",  my_names)
  
  #matrix operation
  my.mx= lapply(mylist, function(x) {if(any(class(x)=="data.frame")) as.matrix(x) else x})
  my.out=Reduce('+', my.mx)/length(input_sim) #index calculation
  
  #rasterize
  my.out_list=list()
  for(i in 1:ncol(my.out)){
    rate.sp=as.data.frame(cbind(cc, data.frame(my.out[,i])))
    names(rate.sp)=c("X", "Y", "ProbIntro")
    coordinates(rate.sp) <- ~X + Y
    gridded(rate.sp) <- TRUE
    rate.sp=raster(rate.sp)
    my.out_list[[i]]=rate.sp
  }
  my.out_list=stack(my.out_list)
  names(my.out_list)=my_names
  return(my.out_list)
}

## Derive abundance CI for each life stage and in each day
indiv_abund_ci <- function(outl=NA,st=1,cores=1,days=0,  breaks=NULL){
  #breaks= c(0.25,0.50,0.75) or c(0.025,0.50,0.975)
  out <- apply(do.call(rbind.data.frame,mclapply(outl, function(x) {
    lapply(1:days, function(y) {
      if(y<=length(x)) {sum(x[[y]][st,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=breaks,na.rm=T); #probs=c(0.025,0.50,0.975)
  colnames(out)<-NULL
  outo <- rbind.data.frame(out,
                           stage=rep(st,nrow(out)),
                           day=as.factor(1:ncol(out)))
  return(t(outo))
}

## SPATIAL Derive abundance CI for each life stage and in each day
indiv_abund_ci_spatial = function(coords=NULL, input_sim=NULL, eval.date=NULL, st=1,  breaks=NULL){
  
  #subset life stage
  mylist=lapply(simout, function(x){Map(function(z, y) z[y, ], x, st)})
  
  #bind together all the days
  mylist=lapply(mylist, function(x){
    pippo=lapply(x, function(y) {data.frame(y)});
    pippo=do.call(cbind, pippo);
  })
  
  #add days name
  col_names=paste0("d_", 1:ncol(mylist[[1]]))
  mylist=lapply(mylist, setNames, col_names)
  
  #subset days
  my_names=col_names[eval.date]
  mylist = lapply(mylist, "[",  my_names)
  
  #matrix operation
  my.mx= lapply(mylist, function(x) {if(any(class(x)=="data.frame")) as.matrix(x) else x})
  my.out <- apply(simplify2array(my.mx), 1:2, quantile, prob =breaks) #first dimension are the breaks; second the grid cells, third dimensions the day
  
  #rasterize
  my.out_list=list()
  for(i in 1:length(eval.date)){
    # i=1
    # names(coords)=c("X", "Y")
    rate.sp=data.frame(t(my.out[, , i]))
    my_quantile_rast=list()
    
    for(j in 1:length(breaks)){ #qui cè l'errore
      tmp=rate.sp[, j]
      tmp=as.data.frame(cbind(cc,tmp))
      names(tmp)[1:2]=c("X", "Y")
      coordinates(tmp) <- ~X+Y
      gridded(tmp) <- TRUE
      tmp=raster(tmp)
      names(tmp)=as.character(breaks[j])
      my_quantile_rast[[j]]=tmp
    }
    my.out_list[[i]]=stack(my_quantile_rast)
  }
  #the output is a list of rasterstack
  names(my.out_list)=my_names
  return(my.out_list)
}

# Derive median and interquartile number of invaded cells per day
invaded_cells_ci <- function(outl=NA,days=0,breaks=NULL){
  #breaks= c(0.25,0.50,0.75) or c(0.025,0.50,0.975)
  out <- apply(do.call(rbind.data.frame,lapply(outl, function(x) {
                lapply(1:days, function(y) {
      if( y<=length(x) ) {
        length(which(x[[y]]>0))
      } else {NA}})})),2,quantile,probs=breaks,na.rm=T)
  colnames(out)<- NULL
  outo <- rbind.data.frame(out,
                           day=as.factor(1:ncol(out)))
  return(t(outo))
}

#spatial rates of mechanistic traits
get_rates_spatial = function(coords=NULL, daily_temp=NULL, species=NULL, rate_fun=NULL, spatial=FALSE, rate=TRUE){
  
  rate_fun = match.fun(rate_fun)
  rate.v <- apply(daily_temp,2, function(x) {
    rate_fun(x, sp=species)})
  
  if(rate==FALSE){
    rate.v=round(1/rate.v,2)
  }
  
  if(spatial) rate.ras <- stack(apply(rate.v,2, function(x) {
    rate.sp <- cbind.data.frame(X=coords[,1],Y=coords[,2],Var=x)
    coordinates(rate.sp) <- ~X+Y
    gridded(rate.sp) <- TRUE
    raster(rate.sp)
  })) else(rate.v)
}
