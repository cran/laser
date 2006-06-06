\name{gamStat}
\alias{gamStat}

\title{ gamStat }
\description{
  Computes the Gamma Statistic of Pybus and Harvey (2000)
}
\usage{
gamStat(x)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
}

\value{
  a list containing the following elements:
  \item{gamstat }{the calculated gamma statistic}
  \item{pval }{One-tailed p-value }
 
}
\references{ Pybus, O. G., Harvey, P. H. 2000.  New uses for incomplete phylogenies. }
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\examples{
  data("plethodon")
  pleth <- plethodon
  gamStat(pleth)
}
\keyword{ univar }
\keyword{ htest }
