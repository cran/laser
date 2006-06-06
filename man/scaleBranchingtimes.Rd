\name{scaleBranchingtimes}
\alias{scaleBranchingtimes}

\title{ Scale Branching Times }
\description{
  Scales a vector of branching times by a known basal divergence time.
}
\usage{
scaleBranchingtimes(x, basal = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{basal}{ estimated time of the basal divergence for the clade }
}
\details{
  This function can be used when analyzing branching times for which it is possible
  to estimate the basal divergence time.  If branch lengths in a specified phylogenetic tree
  are given in units of genetic distance, but the basal divergence is estimated at Z mya, call 
  \samp{res <- scaleBranchingtimes(x, basal = Z)}.  Rate parameters estimated by fitting variants of birth-death
  models will then be in units of lineage births or deaths per million years. 
}
\value{
  a vector of scaled branching times
}

\author{ Dan Rabosky \email{ DLR32@cornell.edu } }

\examples{
  data(plethodon)
  # suppose plethodon basal divergence occurred 23 mya:
  svec <- scaleBranchingtimes(plethodon, basal = 23)
  plotLtt(svec)
  # plots lineages through time in units of time, rather than genetic
  # distance
  
  pureBirth(svec)
  # returns speciation rate in units of lineages/million years
}
\keyword{ manip }
