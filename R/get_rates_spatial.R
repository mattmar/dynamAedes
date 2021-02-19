get_rates_spatial = function(coords=NULL, daily_temp=NULL, species=NULL, rate_fun=NULL, spatial=FALSE, rate=TRUE, np=1){

	rate_fun = match.fun(rate_fun)
	rate.v <- mclapply(apply(daily_temp,2,list), function(x) {
		rate_fun(x, sp=species)},mc.cores=np)

	if(rate==FALSE){
		rate.v=lapply(rate.v, function(x) round(1/x))
	}

	if(spatial) {
		rate.ras <- stack(mclapply(rate.v,function(x) {
			rate.sp <- cbind.data.frame(X=coords[,1],Y=coords[,2],Var=x)
			sp::coordinates(rate.sp) <- ~X+Y
			sp::gridded(rate.sp) <- TRUE
			raster(rate.sp)
		}, mc.cores=np))
		names(rate.ras) <- paste0("d",1:dim(rate.ras)[3])
		return(rate.ras)
	} else(return(do.call(cbind.data.frame, rate.v)))
}