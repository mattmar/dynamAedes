psi_sp = function(coords=NULL, input_sim=NULL, eval_date=NULL, np=1) {
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes")
	} else {
		if( max(eval_date) > max(sapply(input_sim,length)) ) {
			stop("eval_date > than maximum number of simulated days across all iterations...")
		} else {
			mylist <- mclapply(input_sim, function(x){
				mydf <- lapply(x, function(y) {data.frame(y); data.frame("tot_individuals"= colSums(y))});
				mydf <- do.call(cbind, mydf)
			}, mc.cores=np)
			maxl <- max(sapply(mylist,ncol))
			mylist <- mclapply(mylist, function(x) {
				if(length(x)<maxl) {
					cbind.data.frame( x, cbind.data.frame( matrix(rep(rep(0,nrow(x)),(maxl-length(x))), ncol=(maxl-length(x)))) )
				} else {
					x
				}
			},mc.cores=np)
  			#add days and convert to 1 all pixel having a number of individuals >1
			col_names <- paste0( "d_", 1:maxl )
			mylist <- lapply(mylist, setNames, col_names)
			mylist <- lapply(mylist, function(x) {x[x>0] <- 1; return(x)})
  			#subset
			my_names=col_names[eval_date]
			mylist = lapply(mylist, "[",  my_names)
  			#matrix operation
			my.out=Reduce('+', mylist)/length(input_sim)
  		#rasterize
		}
		return(stack(apply(my.out, 2, function(x) {
			rate.sp=as.data.frame(cbind(coords, data.frame(x)))
			names(rate.sp)=c("X", "Y", "ProbIntro")
			rate.sp<-rasterFromXYZ(rate.sp)
			return(rate.sp)
		})))
	}
}
