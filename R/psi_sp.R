#' Probability of successful introduction (spatial)
#'
#' Compute the proportion of successful introductions per each cell of the grid.
#' @param input_sim matrix. \code{dynamAedes.m} compressed output matrix (\code{compressed=TRUE}).
#' @param eval_date positive integer. Define the day(s) to calculate the proportion of successful introductions which should match the column number of the temperature matrix used to inform the model.
#' @param n.clusters positive integer. Define the number of parallel processes.
#' @return \code{psi_sp} returns a raster with the proportion of model iterations that resulted in a viable mosquito population at a given date for a given life stage in each cell of the grid.
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

psi_sp = function(input_sim=NULL, eval_date=NULL, n.clusters=1) {
	if(!is.numeric(eval_date)) {
		stop("eval_date not defined, exiting...")
	}	

	coords = input_sim@coordinates
	input_sim = input_sim@simulation
	
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes")
	} else {
		if( max(eval_date) > max(sapply(input_sim,length)) ) {
			stop("eval_date > than maximum number of simulated days across all iterations...")
		} else {
			mylist <- mclapply(input_sim, function(x){
				mydf <- lapply(x, function(y) {data.frame(y); data.frame("tot_individuals"= colSums(y))});
				mydf <- do.call(cbind, mydf)
			}, mc.cores=n.clusters)
			maxl <- max(sapply(mylist,ncol))
			mylist <- mclapply(mylist, function(x) {
				if(length(x)<maxl) {
					cbind.data.frame( x, cbind.data.frame( matrix(rep(rep(0,nrow(x)),(maxl-length(x))), ncol=(maxl-length(x)))) )
				} else {
					x
				}
			},mc.cores=n.clusters)
  			#add days and convert to 1 all pixel having a number of individuals >1
			col_names <- paste0( "d_", 1:maxl )
			mylist <- lapply(mylist, setNames, col_names)
			mylist <- lapply(mylist, function(x) {x[x>0] <- 1; return(x)})
  			#subset
			my_names=col_names[eval_date]
			mylist = lapply(mylist, "[",  my_names)
  			#matrix operation
			my.out=Reduce('+', mylist)/length(input_sim)
		}
		rateList <- apply(my.out, 2, function(x) {
			rate.sp=as.data.frame(cbind(coords, data.frame(x)))
			#names(rate.sp)=c("X", "Y", "ProbIntro")
			rate.sp<-terra::rast(rate.sp, type="xyz")
			#return(rate.sp)
		})
		rateList <-do.call(c, unname(rateList))
		names(rateList) <- my_names
		return(rateList)
	}
}
