#' Number of invaded cells
#'
#' Compute a summary of the number of invaded cells over model iterations
#' @param input_sim matrix. \code{dynamAedes.m} compressed output matrix (\code{compressed=TRUE}).
#' @param eval_date numeric. Define the day of evaluation; it refers to the column number of the input temperature matrix. 
#' @param breaks numeric vector. Quantile breaks, default interquartile range: \code{c(0.25,0.5,0.75)}.
#' @return \code{icci} returns quantiles of the distribution of the invaded cell number for the specified. The output should be interpreted according to model spatial scale (i.e. \code{scale='rg'} or \code{scale='lc'} give different interpretation).
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

icci <- function(input_sim=NA, eval_date=NULL, breaks=c(0.25,0.5,0.75)){
	if(!is.numeric(eval_date)) {
		stop("eval_date not defined, exiting...")
	}	


	input_sim = input_sim@simulation
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > number of simulated days...")
	}
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes.m()")
	} else {
		out <- apply(do.call(rbind.data.frame,lapply(input_sim, function(x) {
			lapply(eval_date, function(y) {
				if( y<=length(x) ) {
					length(which(colSums(x[[y]])>0))
				} else {NA}})})),2,quantile,probs=breaks,na.rm=T)
		colnames(out)<- NULL
		outo <- rbind.data.frame(out,
			day=as.factor(1:ncol(out)))
	}
	return(t(outo))
}