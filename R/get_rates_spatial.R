get_rates_spatial = function(coords=NULL, daily_temp=NULL, species=NULL, rate_fun=NULL, spatial=FALSE, rate=TRUE, np=1){

	rate_fun = match.fun(rate_fun)
	rate.v <- mclapply(apply(daily_temp,2,list), function(x) {
		rate_fun(x, sp=species)},mc.cores=np)

	if(rate==FALSE){
		rate.v=lapply(rate.v, function(x) round(1/x))
	}

	if(spatial) rate.ras <- stack(mclapply(rate.v,function(x) {
		rate.sp <- cbind.data.frame(X=coords[,1],Y=coords[,2],Var=x)
		coordinates(rate.sp) <- ~X+Y
		gridded(rate.sp) <- TRUE
		raster(rate.sp)
		return(rate.ras)
	}, mc.cores=np)) else(return(do.call(cbind.data.frame, rate.v)))
}