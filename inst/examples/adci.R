\dontrun{
# Make a toy temperature time series
w <- matrix(seq(20,25,length.out=5),ncol=5)*1000

# Run the model
mod <- dynamAedes.m(
	species="koreicus", 
	scale="ws",
	intro.eggs=10, 
	ihwv=2, 
	temps.matrix=w,
	startd=2, 
	endd=5,
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE)

# Derive interquarile of population egg abundance per cell
adci(mod, eval_date=2:4, breaks=c(0.25,0.50,0.75), st=1)
}