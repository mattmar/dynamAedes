\name{psi}
\alias{psi}
\title{
Probability of a successful introduction}
\description{
Derive the proportion of successful introductions over the number of iterations (probability of introduction) for a certain date or range of dates.}
\usage{
psi(input_sim=NULL, eval_date=NULL, stage=0)
}
\arguments{
\item{input_sim}{input data.}
\item{eval_date}{input data.}
\item{stage}{input data.}
}
\value{
When \emp{stage=0} a population is considered exinct when there is no individual (in any stage) in any of the cells of the arena on the evaluated date. If stage!=0 then extinction is relative to the stage compartment. 
}
\author{
Daniele Da Re \email{daniele.dare@uclouvain.be} \cr
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
}
\keyword{internal}
