\name{mccrTest.Rd}
\alias{mccrTest}

\title{ Conduct monte carlo constant rates (gamma statistic) test }
\description{
  Conducts gamma statistic analysis (Pybus and Harvey 2000) for incompletely sampled phylogenies. Phylogenies 
  are simulated to the full clade size under the null hypothesis (constant rate pure birth diversification process) and 
  taxa are randomly pruned from the tree to mimic incomplete sampling. The null distribution of the gamma statistic
  is then tabulated from these phylogenies.
}    
\usage{
mccrTest(CladeSize, NumberMissing, NumberOfReps, ObservedGamma = NULL, fname=NULL)

}


\arguments{
  \item{CladeSize}{  The TRUE clade diversity}
  \item{NumberMissing}{ The number of missing species, e.g., \code{CladeSize} minus the number of taxa in your tree}
  \item{NumberOfReps}{ The number of Monte Carlo simulations to conduct. Recommend >= 5000 }
  \item{ObservedGamma}{ The observed gamma statistic value for the empirical tree. Optional. }
  \item{fname}{ An optional filename where simulated trees are stored, if you generated the trees in another program, such as phylogen}
}
\details{
   If \code{ObservedGamma} is supplied, \code{mccrTest} returns the p-value. 
}
\value{
   A list with the following components:
  \item{null.gamma}{ The null distribution of gamma. You can plot a histogram or otherwise inspect these values... }
  \item{critical.value}{ The 0.05 percentile of the null distribution. This is the value corresponding to alpha = 0.05 }
  \item{p.value}{ The actual p-value, only returned if \code{ObservedGamma} is supplied by user}
}
\references{ 

Pybus, O. G. and Harvey, P. H. 2000. Testing macro-evolutionary models using incomplete molecular phylogenies. Proceedings of the Royal Society of London. Series B. Biological Sciences, 267, 2267-2272.

}
\author{ Dan Rabosky \email{drabosky@umich.edu} }

\seealso{ \code{\link{gamStat}} }
\examples{
 	mccrTest(CladeSize=25, NumberMissing=5, NumberOfReps=50);
}
\keyword{ univar }
\keyword{ htest }




