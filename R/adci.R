#' Summaries of mosquito abundance
#'
#' Summaries of mosquito abundance at each life stage for each day 
#' @param input_sim (matrix) dynamAedes compressed output matrix 
#' @param eval_date (integer) define the day of successful introduction evaluation, referring to the column number of the temperature matrix used to inform the model. 
#' @param stage (integer) 0 (all), 1 (egg), 2 (juvenile), 3 (adult), 4 (diapausing egg).
#' @param n.clusters (numeric) define the number of parallel processes.
#' @param breaks quantile breaks, default the first, the second and the third c(0.25,0.5,0.75)
#' @return Returns a table with the summary of mosquito abundance at each life stage for each day. 
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

adci <- function(input_sim=NA, stage=1, n.clusters=1, eval_date=0, breaks=c(0.25,0.5,0.75)){
	if(stage<1|stage>4) {
		stop("stage can be 1 (egg), 2 (juvenile), 3 (adult), 4 (diapausing eggs) ...")
	}
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > than number of simulated days...")
	} else {
		out <- apply(do.call(rbind.data.frame,mclapply(input_sim, function(x) {
			lapply(eval_date, function(y) {
				if(y<=length(x)) {sum(x[[y]][stage,],na.rm=T)} else {NA}})},mc.cores=n.clusters)),2,quantile,probs=breaks,na.rm=T)
		colnames(out)<-NULL
		outo <- rbind.data.frame(out,
			stage=rep(stage,nrow(out)),
			day=as.factor(1:ncol(out)))
		return(t(outo))
	}
}