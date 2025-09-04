\donttest{
# Make toy temperature time series and coordinates
w <- matrix(as.integer(seq(25,35,length.out=250)*1000),ncol=10)
cc <- matrix(c(seq(9.5,10,length.out=25),seq(41,42,length.out=25)),ncol=2)

# Run the model
mod <- dynamAedes.m(
	species="albopictus", 
	scale="rg",
	intro.eggs=100, 
	jhwv=10, 
	temps.matrix=w,
	cells.coords=cc,
	coords.proj4="4326",  
	startd="2021-06-21", 
	endd="2021-06-25",
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE,
	seeding=TRUE)

# Interquarile of population egg abundance per cell
adci_sp(mod, coords=cc, eval_date=3, breaks=c(0.025,0.5,0.975), st=3)
}