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

