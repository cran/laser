\name{plotLtt}
\alias{plotLtt}

\title{ plotLtt }
\description{
  Plots log-lineages through time given a vector of branching times.
}
\usage{
plotLtt(x)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
}

\value{
  a plotting function with no return value
}

\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\examples{

data("plethodon")
pleth.btimes <- plethodon
plotLtt(pleth.btimes)

}
\keyword{ aplot }
\keyword{ hplot }
