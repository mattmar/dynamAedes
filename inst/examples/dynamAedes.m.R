\donttest{
## Run dynamAedes at local scale for 5 days
# Make a toy temperature time series
w <- matrix(seq(20,25,length.out=5),ncol=5)*1000
# Run the model
	simout <- dynamAedes.m(
	species="koreicus", 
	scale="ws",
	intro.eggs=10, 
	jhwv=2, 
	temps.matrix=w, 
	startd="2021-06-21", 
	endd="2021-06-25",
	lat=42,
	long=8,
	n.clusters=1, 
	iter=1,
	compressed.output=TRUE,
	verbose=FALSE)
}