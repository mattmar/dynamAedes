\donttest{
# Make toy temperature time series and coordinates
w <- matrix(as.integer(seq(25,35,length.out=250)*1000),ncol=10)
cc <- matrix(c(rep(seq(1,5),5),rev(rep(seq(1,5),each=5))),ncol=2)

# Run the model
mod <- dynamAedes.m(
	species="albopictus", 
	scale="rg",
	intro.eggs=10, 
	jhwv=10, 
	temps.matrix=w,
	cells.coords=cc,  
	startd="2021-06-21", 
	endd="2021-06-25",
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE,
	seeding=TRUE)

# Proportion of viable populations on the last day of simulations per cell
psi_sp(coords = cc, input_sim = mod, eval_date = 3)
}