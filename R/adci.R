#' Summaries of mosquito abundance.
#'
#' Summaries of mosquito abundance at each life stage for each day.
#' @param input_sim matrix. dynamAedes.m compressed output matrix.
#' @param eval_date positive integer. Define the day of successful introduction evaluation, referring to the column number of the temperature matrix used to inform the model. 
#' @param stage positive integer. 0 (all), 1 (egg), 2 (juvenile), 3 (adult), 4 (diapausing egg).
#' @param n.clusters positive integer. Define the number of parallel processes.
#' @param breaks numeric vector. Quantile breaks, default the first, the second and the third quantile: \code{c(0.25,0.5,0.75)}.
#' @return Returns a table with the summary of mosquito abundance per life stage (or substage if \code{compressed.output=FALSE} in \code{dynamAedes.m} function) for each day. 
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

adci <- function(input_sim=NA, stage=1, n.clusters=1, eval_date=0, breaks=c(0.25,0.5,0.75)){
	if(attributes(input_sim)$compressed) {
		if(stage<1|stage>4) {
			stop("stage can be 1 (egg), 2 (juvenile), 3 (adult), 4 (diapausing eggs) ...")
		}
		if( max(eval_date) > max(sapply(input_sim,length)) ) {
			stop("eval_date > than number of simulated days...")
		} else {
			out <- apply(do.call(rbind.data.frame,mclapply(input_sim, function(x) {
				lapply(eval_date, function(y) {
					if(y<=length(x)) {sum(as.data.frame(x[[y]])[stage,],na.rm=T)} else {NA}})},mc.cores=n.clusters)),2,quantile,probs=breaks,na.rm=T)
			colnames(out)<-NULL
			outo <- rbind.data.frame(out,
				stage=rep(stage,nrow(out)),
				day=as.factor(1:ncol(out)))
			return(t(outo))
		}
	} else {
		if( max(eval_date) > max(sapply(input_sim,length)) ) {
			stop("stage can be 1 (egg), 2 (juvenile), 3 (adult), 4 (diapausing eggs) ...")
		} else {
			tmp <- mclapply(1:length(input_sim), function(x) {as.matrix(lapply(eval_date, function(dt) {input_sim[x][[1]][dt][[1]][stage,,]}))},mc.cores=n.clusters)
			tmp1 <- lapply(eval_date, function(dt) cbind.data.frame(
				apply(simplify2array(lapply(tmp,function(x) (x[[dt]]))),1,quantile,breaks),
				day=dt,
				stage=stage))
			outo <- do.call(rbind.data.frame,lapply(tmp1, function(x) {x$q <- row.names(x);row.names(x)<-NULL;return(x)}))
			return(outo)
		}
	}
}