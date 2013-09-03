\name{gamStat}
\alias{gamStat}

\title{ gamStat }
\description{
  Computes the Gamma Statistic of Pybus and Harvey (2000), using branching times
}
\usage{
gamStat(x, return.list=TRUE)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{return.list}{ logical: return a list with pvalue and the statistic? }
}

\value{
  If \code{return.list} = TRUE, a list containing the following elements:
  \item{gamstat }{the calculated gamma statistic}
  \item{pval }{One-tailed p-value }
 
 If \code{return.list} = FALSE, it simply returns a single number, the gamma statistic 
}
\references{ Pybus, O. G., Harvey, P. H. 2000.  New uses for incomplete phylogenies. }
\author{ Dan Rabosky \email{drabosky@umich.edu} }

\examples{
  data("plethodon")
  pleth <- plethodon
  gamStat(pleth)
}
\keyword{ univar }
\keyword{ htest }
