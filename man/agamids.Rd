\name{agamids}
\alias{agamids}
\docType{data}
\title{ Ultrametric Phylogeny of Australian Agamid Lizards}
\description{
  Ultrametric phylogenetic tree for 69 Australian Agamid lizards.  Tree was constructed 
  by maximum likelihood from 1800 bp of mtDNA using the GTR + I + G model of sequence evolution and made
  ultrametric using Non-Parametric Rate Smoothing.
}
\usage{data(agamids)}

\source{
  Harmon, L. J., J. A. Schulte, A. Larson, and J. B. Losos.  2003.  
  Tempo and mode of evolutionary radiation in Iguanian lizards. \emph{Science} 301:961-964.
}
\references{
Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 
}
\examples{
data(agamids)
agbtimes <- getBtimes(string = agamids)
gamStat(agbtimes)
plotLtt(agbtimes)

}
\keyword{datasets}
