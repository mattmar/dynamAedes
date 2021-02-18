psi <- function(input_sim=NULL, eval_date=NULL, stage=0){
	if(stage<0|stage>3) {
		stop("stage can be 0 (all), 1 (egg), 2 (juvenile), 3 (adult) ...")
	}
	st <- c("all","egg","juvenile","adult","diapausing egg")
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > than maximum number of simulated days across all iterations...")
	} else {
		pe_out <- sapply(eval_date, function(y) {
			pe <- sum(unlist(lapply(input_sim, function(x) {
				if(stage==0) {
					length(which(sum(unlist(x[y]))>0));
				} else {
					length(which(x[y][[1]][stage,]>0))
				}
			}))) / length(input_sim)
			return(pe)
		})

		pe_out <- data.frame(unlist(pe_out))
		names(pe_out) <- "p_success"
		rowid <- paste0("Day ", eval_date)
		pe_out <- cbind.data.frame("Days_after_intro"=rowid, "p_success"=pe_out)
		pe_out$stage <- rep(st[stage+1], nrow(pe_out))
	}
	return(pe_out)
}