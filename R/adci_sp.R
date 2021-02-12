adci_sp = function(coords=NULL, input_sim=NULL, eval.date=NULL, st=1,  breaks=c(0.25,0.5,0.75)){
	if( all(unlist(lapply(input_sim, function(x) { sapply(x,length) } ))==4) ) {
		stop("Non-spatial data, set scale='lc' or scale='rg' in dynamAedes")
	} else {
    	# subset life stage
		mylist=lapply(input_sim, function(x){ Map(function(z, y) z[y, ], x, st) })
    	# bind together all the days
		mylist=lapply(mylist, function(x){
			pippo=lapply(x, function(y) {data.frame(y)});
			pippo=do.call(cbind, pippo);
		})
    	# add days name
		col_names=paste0("d_", 1:ncol(mylist[[1]]))
		mylist=lapply(mylist, setNames, col_names)
    	# subset days
		my_names <- col_names[eval.date]
		mylist <- lapply(mylist, "[",  my_names)
    	# matrix operation
		my.mx <- lapply(mylist, function(x) {if(any(class(x)=="data.frame")) as.matrix(x) else x})
		my.out <- apply(simplify2array(my.mx), 1:2, quantile, prob=breaks)
    	# rasterize
		outr <- lapply(1:length(eval.date), function(x) {
			rate.sp <- data.frame(t(my.out[, , x]))

			rb <- lapply(1:length(breaks),function(y) {
				tmp=rate.sp[, y]
				tmp=as.data.frame(cbind(cc,tmp))
				names(tmp)[1:2]=c("X", "Y")
				coordinates(tmp) <- ~X+Y
				gridded(tmp) <- TRUE
				tmp=raster(tmp)
				names(tmp)=as.character(breaks[y])
				return(tmp)
			})
			return(stack(rb))
		})
	}
	return(outr)
}