\name{yuleWindow}
\alias{yuleWindow}
\title{ Fit Yule Model (Pure Birth) to Temporal Window of Branching Times }
\description{
  This function fits the Yule model to any temporal window that includes at least 
  one branching time.  
}
\usage{
yuleWindow(x, st1 = x[1], st2 = 0)
}
\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{st1}{ the start of the interval you wish to examine }
  \item{st2}{ the end of the interval you wish to examine }
}
\details{
  'st1' and 'st2' are given in divergence units before present.  If we had a set 
  of branching times with an initial divergence 100 mya, \samp{yuleWindow(x, 75, 25)} would
  fit the pure birth model to the portion of the tree between 75 and 25 mya.  
  Calling \samp{yuleWindow(x, x[1], 0)} will fit the model to the entire tree (identical 
  to \samp{pureBirth(x)}).
  
  Note that st1 must be greater than st2, because they are given in units of divergence 
  before present.  
}
\value{
  a list with the following components:
  \item{LH }{the log-likelihood at the maximum}
  \item{smax }{the speciation rate giving the maximum log-likelihood}
 
}
\references{ 
Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed evolutionary 
process. \emph{Philos. Trans. R. Soc. Lond. B} 344:305-311.

Nee, S.  2001.  Inferring speciation rates from phylogenies. 
\emph{Evolution} 55:661-668.

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{
  This can be used in conjunction with other models in this package to 
  test \emph{a priori} hypotheses of rate variation.  If, for example, it was
  hypothesized that a particular climatic event occurring T mya shifted diversification rates,
  you can use \samp{yuleWindow(x, x[1], T)} and \samp{yuleWindow(x, T, 0)} to obtain
  the log-likelihoods for these two temporal windows.  The log-likelihoods can then be summed
  to obtain the likelihood of the full set of branching times under a 2-rate Yule model with 
  an \emph{a priori} hypothesized rate-shift.  This model would only have two free parameters, in contrast
  to the \code{\link{yule2rate}} model, with three parameters (and thus, the AIC would be computed 
  as (- 2 * (sum of log-likelihoods) + 4)). 

}
\seealso{ \code{\link{fitdAICrc}}, \code{\link{pureBirth}}, \code{\link{rvbd}}, \code{\link{bd}} }
\examples{
  data(agamids)
  agbtimes <- getBtimes(string = agamids)
  yuleWindow(agbtimes, 0.22, 0.10)
  # fits Yule model to temporal window between 0.22 and 0.10 divergence
  # units before present
  
}
\keyword{ models }
\keyword{ htest }
