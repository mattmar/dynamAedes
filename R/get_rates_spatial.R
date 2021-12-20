#' Proportion of successful introductions
#'
#' Compute the proportion of "successful" introductions.
#' @param coords matrix (optional). A matrix reporting the spatial coordinates of the temperature observations.
#' @param temps.matrix matrix, vector. A matrix or a vector of daily average temperature (in Celsius degree /1000) used to fit the life cycle rates. This matrix must be organised with the daily temperature observations as columns and the geographic position of the i-grid cell as rows.
#' @param species character. Select what species to model: \code{"aegypti"}, \code{"albopictus"}, \code{"japonicus"}, \code{"koreicus"}. Default \code{species = "aegypti"}. 
#' @param rate_fun character. Define one of the temperature-dependent functions of the model, e.g. ".e.hatch.rate". 
#' @param spatial logical. Get a raster as output. Default \code{FALSE}
#' @param rate logical. Returns the daily rate. Default \code{TRUE}
#' @param n.clusters postive integer. Define the number of parallel processes.
#' @return Returns the estimated value of a given temperature-dependent function. If \code{spatial = TRUE}, it will return a raster. 
#' @author Matteo Marcantonio \email{marcantoniomatteo@gmail.com}, Daniele Da Re \email{daniele.dare@uclouvain.be}
#' @export

get_rates_spatial <- function(coords=NULL, temps.matrix=NULL, species="aegypti", rate_fun=NULL, spatial=FALSE, rate=TRUE, n.clusters=1){

  rate_fun = match.fun(rate_fun)

	if (ncol(as.matrix(temps.matrix))==1) {
	  inTemp=list(temps.matrix)
	  rate.v <- mclapply(inTemp, function(x) { rate_fun(unlist(x)/1000, sp=species)},mc.cores=n.clusters)  
	}
	
	if (ncol(as.matrix(temps.matrix))>=2) {
	  inTemp=apply(temps.matrix,2,list)
	  rate.v <- mclapply(inTemp, function(x) { rate_fun(unlist(x)/1000, sp=species)},mc.cores=n.clusters)
	} 

	if (rate==FALSE) {
		rate.v=lapply(rate.v, function(x) round(1/x))
	}

	if(spatial) {
		rate.ras <- stack(mclapply(rate.v, function(x) {
			rate.sp <- cbind.data.frame(X=coords[,1],Y=coords[,2],Var=x)
			rate.sp<-rasterFromXYZ(rate.sp)
		}, mc.cores=n.clusters))
		names(rate.ras) <- paste0("d",names(as.data.frame(temps.matrix)))
		return(rate.ras)
	} else(return(do.call(cbind.data.frame, rate.v)))
}
