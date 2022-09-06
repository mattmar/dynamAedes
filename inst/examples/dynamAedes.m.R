## Run dynmAedes at local scale for 5 days
# Make a toy temperature time series
\dontrun{
w <- matrix(seq(20,25,length.out=5),ncol=5)*1000
# Run the model
	dynamAedes.m(
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
}