\name{gamStat.batch}
\alias{gamStat.batch}
\alias{batchgam}

\title{ Compute Gamma Statistic For Batch of Branching Times }
\description{
  Takes a matrix or dataframe of branching times obtained from a 
  set of phylogenies using \code{\link{getBtimes.batch}} and performs the 
  MCCR test of Pybus and Harvey (2000).  Can also be used to generate the posterior
  distribution of the gamma statistic for Bayesian inference.

}
\usage{
gamStat.batch(x, stat = NULL)
}

\arguments{
  \item{x}{ a matrix or dataframe containing branching times for phylogenies }
  \item{stat}{ value of the gamma statistic for the test phylogeny }
}
\details{
  'x' is obtained through use of \code{\link{getBtimes.batch}}.  
  
  'stat' is the calculated value of the gamma statistic you wish to compare
  to the null distribution for a given level of sampling and number of tips.  If 
  \samp{stat = NULL} (the default), the function returns a list with calculated 
  gamma statistics for each set of branching times, as well as the critical value of the 
  distribution (alpha = 0.05).  Observed gamma statistic values less than the critical value
  have probability < 0.05 under the null hypothesis that rates have not decreased over time.  If 
  \samp{stat} is specified, the function will determine the p-value of the observed value under 
  the null hypothesis by tabulating the number of simulated values above and below the observed 
  value. 
}
\value{
  a list with the following components:
  \item{gamstat }{ a vector of calculated gamma statistics for the matrix of branching times}
  \item{gamcrit }{ the critical value of the distribution specified by the simulated trees (alpha = 0.05)}
  \item{pval}{ the Type I error probability of the observed value under the null distribution 
    specified by the input data \samp{x}}
}
\references{ 
Pybus, O. G., and P. H. Harvey.  2000.  Testing macro-evolutionary models using 
incomplete molecular phylogenies. \emph{Proc. R. Soc. Lond. Biol. B} 267:2267-2272.
}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ 
  This performs the MCCR test of Pybus and Harvey (2000).  However, you must
  generate phylogenies under a rate-constant pure-birth (Yule) model using either
  \emph{Phyl-O-Gen} (\url{http://evolve.zoo.ox.ac.uk/software/PhyloGen/main.html}) or another 
  program.  
  
  This function can also be used for Bayesian inference.  For example, you can tabulate the distribution 
  of the gamma statistic from the posterior distribution of phylogenetic trees 
  sampled using MCMC or MCMCMC.  
  
  
}
\seealso{ \code{\link{gamStat}}, \code{\link{getBtimes.batch}}, \code{\link{fitdAICrc.batch}} }
\examples{

data("rtrees50")
  write.table(rtrees50, file = 'temp.txt', quote=FALSE, 
                row.names = FALSE, col.names = FALSE)
#generates file 'temp.txt' with data formatted appropriately 
intrees <- getBtimes.batch(file = 'temp.txt')
result <- gamStat.batch(intrees)
hist(result$gamstat)
result$gamcrit
unlink('temp.txt') #delete temp file

}
\keyword{ univar }
\keyword{ htest }
