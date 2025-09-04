#' Proportion of successful introductions
#'
#' Compute the proportion of "successful" introductions.
#' @param input_sim matrix. \code{dynamAedes.m} compressed output matrix (\code{compressed=TRUE}). 
#' @param eval_date positive integer. define the day(s) to calculate the proportion of successful introductions which should match the column number of the temperature matrix used to inform the model.
#' @return \code{psi} returns the proportion of model iterations that resulted in a viable mosquito population (defined as: iterations with at least one individual alive in any life stage) at a given date.
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

psi <- function(input_sim=NULL, eval_date=NULL){
	if(!is.numeric(eval_date)) {
		stop("eval_date not defined, exiting...")
	}

	input_sim = input_sim@simulation
	if( max(eval_date) > max(sapply(input_sim,length)) ) {
		stop("eval_date > maximum number of simulated days across all iterations...")
	} else {
		pe_out <- sapply(eval_date, function(y) {
			pe <- sum(unlist(lapply(input_sim, function(x) {
				length(which(sum(unlist(x[y]))>0));
			}))) / length(input_sim)
			return(pe)
		})

		pe_out <- data.frame(unlist(pe_out))
		names(pe_out) <- "p_success"
		rowid <- paste0("Day ", eval_date)
		pe_out <- cbind.data.frame("Days_after_intro"=rowid, "p_success"=pe_out)
		pe_out$stage <- rep("Population", nrow(pe_out))
	}
	return(pe_out)
}