# Make toy temperature time series and coordinates
\dontrun{
w <- matrix(as.integer(seq(25,35,length.out=250)*1000),ncol=10)
cc <- matrix(c(rep(seq(1,5),5),rev(rep(seq(1,5),each=5))),ncol=2)

# Run the model
mod <- dynamAedes(
	species="albopictus", 
	scale="rg",
	intro.eggs=10, 
	ihwv=10, 
	temps.matrix=w,
	cells.coords=cc,  
	startd=2, 
	endd=10,
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE,
	seeding=TRUE)

# Plot proportion of viable populations on the last day of simulations per cell
plot(psi_sp(coords = cc, input_sim = mod, eval_date = 9))
}