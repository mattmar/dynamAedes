\name{psi_sp}
\alias{psi_sp}
\title{
Probability of successful introduction}
\description{
Derive the proportion of successful introductions over the number of iterations (probability of introduction)  per each pixel for a certain date.}
\usage{
psi_sp(coords=NULL, input_sim=NULL, eval_date=NULL)
}
\arguments{
\item{coords}{input data.}
\item{input_sim}{input data.}
\item{eval_date}{input data.}
}
\value{
A population is considered exint when on the evaluated date there is no individual (in any stage) in any of the cells of the arena.
}
\author{
Daniele Da Re \email{daniele.dare@uclouvain.be} \cr
Marcantonio Matteo \email{marcantoniomatteo@gmail.com} \cr
}
\keyword{internal}
