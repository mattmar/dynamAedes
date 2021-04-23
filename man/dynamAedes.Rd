\name{dynamAedes}
\alias{dynamAedes}
\title{
Life cycle `Aedes` mosquitoes
}
\description{
Function to simulate population dynamics of `Aedes` mosquitoes}
\usage{
dynamAedes(species="aegypti", intro.eggs=0, intro.adults=0, 
	intro.juveniles=0, scale="ws",intro.cells=NULL, ihwv=1, 
	temps.matrix=NULL, startd=1, endd=10, cells.coords=NULL, 
	lat=0, long=0, road.dist.matrix=NULL, country=NA, 
	intro.year=2020, iter=1, n.clusters=1, cluster.type="PSOCK", 
	sparse.output=FALSE, compressed.output=TRUE, suffix=NA, 
	cellsize=250, maxadisp=600, dispbins=10, verbose=FALSE)
}
\arguments{
\item{species}{input data.}
\item{intro.eggs}{input data.}
\item{intro.adults}{input data.}
\item{intro.juveniles}{input data.}
\item{scale}{input data.}
\item{intro.cells}{input data.}
\item{ihwv}{input data.}
\item{temps.matrix}{input data.}
\item{startd}{input data.}
\item{endd}{input data.}
\item{cells.coords}{input data.}
\item{lat}{input data.}
\item{long}{input data.}
\item{road.dist.matrix}{input data.}
\item{country}{input data.}
\item{intro.year}{input data.}
\item{iter}{input data.}
\item{n.clusters}{input data.}
\item{cluster.type}{input data.}
\item{sparse.output}{input data.}
\item{compressed.output}{input data.}
\item{suffix}{input data.}
\item{verbose}{input data.}
\item{cellsize}{input data.}
\item{maxadisp}{input data.}
\item{dispbins}{input data.}
\item{verbose}{input data.}
}
\value{
Matrix or a list of matrixes containing, for each iteration, the number of individuals in each life stage per day (and for each grid cell of the study area if scale="lc" or "rg". If the argument compressed.output=FALSE (default TRUE), he model returns the daily number of individuals in each life stage sub-compartment.	
}
\author{
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
Daniele Da Re \email{daniele.dare@uclouvain.be} \cr
}
\keyword{internal}
