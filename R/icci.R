#' Number of invaded cells
#'
#' Compute a summary of the number of invaded cells over model iterations
#' @param input_sim (matrix) dynamAedes compressed output matrix 
#' @param eval_date (integer) define the day of successful introduction evaluation, referring to the column number of the temperature matrix used to inform the model. 
#' @param breaks quantile breaks, default the first, the second and the third c(0.25,0.5,0.75)
#' @return icci returns the Interquartile range of the number of invaded cells for the specified day over model all iterations. The output should be interpreted according to dynamAedes spatial scale (i.e. scale='rg' or 'lc').
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

icci <- function(input_sim=NA, eval_date=0, breaks=c(0.25,0.5,0.75)){
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > than number of simulated days...")
	}
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes()")
	} else {
		out <- apply(do.call(rbind.data.frame,lapply(input_sim, function(x) {
			lapply(1:length(eval_date), function(y) {
				if( y<=length(x) ) {
					length(which(colSums(x[[y]])>0))
				} else {NA}})})),2,quantile,probs=breaks,na.rm=T)
		colnames(out)<- NULL
		outo <- rbind.data.frame(out,
			day=as.factor(1:ncol(out)))
	}
	return(t(outo))
}