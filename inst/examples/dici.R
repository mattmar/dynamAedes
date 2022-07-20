\dontrun{
# Make toy temperature time series and coordinates
rr <- raster(vals=NA, nrows=5, ncols=5, xmn=0, xmx=500, ymn=0, ymx=500, resolution=100)
w <- as.data.frame(matrix(as.integer(seq(25,35,length.out=250)*1000),ncol=10))
cc <- coordinates(rr)
dd <- apply(as.matrix(dist(coordinates(cc[3:10,]))),2,function(x) round(x/1000,1)*1000)
colnames(dd) <- 3:10
row.names(dd) <- colnames(dd)
row.names(cc) <- 1:25
row.names(w) <- 1:25


# Run the model (not workng)
mod <- dynamAedes.m(
	species="albopictus", 
	scale="lc",
	intro.juveniles=100, 
	ihwv=10, 
	temps.matrix=w,
	cells.coords=cc,
	road.dist.matrix=dd,
	intro.cells=3,  
	startd=2, 
	endd=8,
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	cellsize=100,
    maxadisp=500,
    dispbins=10,
    seeding=TRUE,
    verbose=TRUE,
  	compressed.output=TRUE
  	)

# Derive estimates of mosquito dispersal of the simulated mosquito populations (only when scale = "lc") for any simulated day (in this case for 9 days from start and end of the simulate period).
dici(mod, coords=cc, eval_date=seq(1,9), breaks=c(0.25,0.50,0.975))
plot(dici(mod, coords=cc, eval_date=seq(1,9), breaks=c(0.25,0.50,0.975), space=TRUE))

}