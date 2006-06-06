\name{bd}
\alias{bd}

\title{ Fit Rate-Constant Birth-Death Model to Branching Times }
\description{
  Finds maximum likelihood estimates of the net diversification rate \emph{r}
  (speciation rate \emph{S} minus the extinction rate \emph{E}) and the extinction
  fraction \emph{a} = \emph{E/S}, using branching times derived from an
  ultrametric phylogenetic tree. 
}
\usage{
bd(x, ai = c(0.1, 0.5, 0.9))
}

\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{ai}{ a vector of initial \emph{a} parameterizations for the optimization algorithm }
}
\details{
  Non-linear optimization can be exceedingly difficult, and the algorithms used here can become trapped
  on local (rather than global) optima.  The default 'ai' parameters specified above fit the constant-rate
  birth-death model to branching times using three initial \emph{a} values.  You should check your results
  against those obtained using the \code{\link{pureBirth}} model.  If the log-likelihood under \samp{bd} is
  less than \samp{pureBirth}, you should explore alternative initial parameterizations.  For example, 
  \samp{ai = seq(0.05, 0.99, length.out = 20)} would attempt the optimization with 20 equally spaced \emph{a} 
  values on the interval (0.05, 0.99).  
  
  I have found the default option to be satisfactory for all phylogenies I have examined.
}
\value{
  a list with the following components:
  \item{LH }{ the log-likelihood at the maximum}
  \item{aic}{ the Akaike Information Criterion }
  \item{r}{ the net diversification rate giving the maximum log-likelihood }
  \item{a}{ the extinction fraction giving the maximum log-likelihood}
}
\references{ 
Kendall, D. G.  1948.  On the generalized "birth-and-death" process. 
\emph{Ann. Math. Stat.} 19:1-15.

Nee, S., E. C. Holmes, R. M. May, and P. H. Harvey.  1994a.  Extinction rates 
can be estimated from molecular phylogenies. \emph{Philos. Trans. R. Soc. Lond. B} 344:77-82.

Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed evolutionary process. 
\emph{Philos. Trans. R. Soc. Lond. B 344:305-311.}

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{pureBirth}}, \code{\link{fitdAICrc}},\code{\link{yule-n-rate}},\code{\link{rvbd}}}
\examples{
  data(plethodon)
  result <- bd(plethodon)  
}
\keyword{ models }
