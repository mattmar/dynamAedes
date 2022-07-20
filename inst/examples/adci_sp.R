\dontrun{
# Make toy temperature time series and coordinates
w <- matrix(as.integer(seq(25,35,length.out=250)*1000),ncol=10)
cc <- matrix(c(rep(seq(1,5),5),rev(rep(seq(1,5),each=5))),ncol=2)

# Run the model
mod <- dynamAedes.m(
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

# Plot interquarile of population egg abundance per cell
plot(adci_sp(mod, coords=cc, eval_date=6, breaks=c(0.025,0.5,0.975), st=3))
}