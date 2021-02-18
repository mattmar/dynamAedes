adci <- function(input_sim=NA, stage=1, cores=1, eval_date=0, breaks=c(0.25,0.5,0.75)){
	if(stage<1|stage>3) {
		stop("stage can be 1 (egg), 2 (juvenile), 3 (adult) ...")
	}
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > than number of simulated days...")
	} else {
		out <- apply(do.call(rbind.data.frame,mclapply(input_sim, function(x) {
			lapply(eval_date, function(y) {
				if(y<=length(x)) {sum(x[[y]][stage,],na.rm=T)} else {NA}})},mc.cores=cores)),2,quantile,probs=breaks,na.rm=T)
		colnames(out)<-NULL
		outo <- rbind.data.frame(out,
			stage=rep(stage,nrow(out)),
			day=as.factor(1:ncol(out)))
		return(t(outo))
	}
}