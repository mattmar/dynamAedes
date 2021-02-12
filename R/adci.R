adci <- function(input_sim=NA, stage=1, cores=1, days=0, breaks=c(0.25,0.5,0.75)){
	out <- apply(do.call(rbind.data.frame,mclapply(input_sim, function(x) {
		lapply(days, function(y) {
			if(y<=length(x)) {sum(x[[y]][stage,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=breaks,na.rm=T); 
	colnames(out)<-NULL
	outo <- rbind.data.frame(out,
		stage=rep(stage,nrow(out)),
		day=as.factor(1:ncol(out)))
	return(t(outo))
}