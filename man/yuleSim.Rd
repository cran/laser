\name{yuleSim}
\alias{yuleSim}

\title{ yuleSim }
\description{
  Simulate Branching Times Under the Pure Birth Model
}
\usage{
yuleSim(ntaxa, nsets, lambda = 0.01)
}

\arguments{
  \item{ntaxa}{ number of taxa in each set of branching times }
  \item{nsets}{ number of datasets of size \code{ntaxa} to be simulated }
  \item{lambda}{ speciation rate (default = 0.01) }
}

\value{
  an n x m matrix of branching times, where n is size \code{ntaxa}
   and m is size \code{ nsets }
}

\author{ Dan Rabosky \email{DLR32@cornell.edu} }


\examples{
  testdata <- yuleSim(25, 50, lambda = .001)
}
\keyword{ datagen }
