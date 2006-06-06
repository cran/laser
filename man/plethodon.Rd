\name{plethodon}
\alias{plethodon}
\docType{data}
\title{ plethodon}
\description{
  Branching times for plethodontid salamanders from Highton and Larson (1979) 
  and analyzed by Nee et al. 1994 and Nee 2001
}
\usage{data(plethodon)}
\format{
  The data are represented by a numeric vector of branching times
}

\source{
  Highton, R., and A. Larson.  1979.  The genetic relationships of the salamanders of the genus
  \emph{Plethodon}.  \emph{Syst. Zool.} 28:579 - 599.
}
\references{
  Nee, S., E. C. Holmes, R. M. May, and P. H. Harvey.  1994a.  
  Extinction rates can be estimated from molecular phylogenies. 
  \emph{Philos. Trans. R. Soc. Lond. B} 344:77-82.
  
  Nee, S.  2001.  Inferring speciation rates from phylogenies. 
  \emph{Evolution} 55:661-668.
}
\examples{
data(plethodon)
plotLtt(plethodon)
}
\keyword{datasets}
