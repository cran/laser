\name{DensityDependent}
\alias{DDX}
\alias{DDL}

\title{ Fit Density Dependent Speciation Model to Branching Times }
\description{
  Functions that fit density-dependent speciation rate models to branching times
  derived from phylogenetic data.  \samp{DDX} and \samp{DDL} fit exponential and logistic 
  variants of the density-dependent speciation rate model.   
}    
\usage{
DDX(x)
DDL(x)
}


\arguments{
  \item{x}{ a numeric vector of branching times }
}
\details{
  \samp{DDX} models the speciation rate as a function of the number of extant 
  lineages at any point in time, r(t) = r0 * (Nt \^ (-x)), where r0 is the initial speciation
  rate, Nt is the number of lineages at some time t, and x is a parameter controlling the magnitude of the 
  rate change.
  
  \samp{DDL} models the speciation rate as r(t) = r0 * (1 - Nt / K), where r0 is the initial speciation
  rate, Nt is the number of extant lineages at some time t, and K is analogous to the 'carrying capacity'
  parameter of population ecology.  
}
\value{
  a list with the following components:
  \item{LH }{ The log-likelihood at the maximum}
  \item{aic }{the Akaike Information Criterion}
  \item{r1 }{ the initial speciation rate}
  \item{kparam }{ the K parameter in the logistic density dependent model}
  \item{xparam }{ the x parameter in the density-dependent exponetial model}
}
\references{ 

Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed evolutionary 
process. \emph{Philos. Trans. R. Soc. Lond. B} 344:305-311.

Nee, S., A. Mooers, and P. H. Harvey.  1992.  Tempo and mode of evolution 
revealed from molecular phylogenies. \emph{Proc. Natl. Acad. Sci. USA} 89:8322-8326.

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{yule-n-rate}}, \code{\link{rvbd}}, \code{\link{fitdAICrc}} }
\examples{
  data(plethodon)
  resX <- DDX(plethodon)
  resL <- DDL(plethodon)
}
\keyword{ models }
