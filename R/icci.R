icci <- function(input_sim=NA, days=0, breaks=c(0.25,0.5,0.75)){
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes()")
	} else {
		out <- apply(do.call(rbind.data.frame,lapply(input_sim, function(x) {
			lapply(1:length(days), function(y) {
				if( y<=length(x) ) {
					length(which(colSums(x[[y]])>0))
				} else {NA}})})),2,quantile,probs=breaks,na.rm=T)
		colnames(out)<- NULL
		outo <- rbind.data.frame(out,
			day=as.factor(1:ncol(out)))
	}
	return(t(outo))
}