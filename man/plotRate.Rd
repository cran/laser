\name{plotRate}
\alias{plotRate}
\title{ Plot speciation and extinction rates through time }
\description{
  Plots diversification rates through time estimated from fits of exponential speciation and extinction models (Rabosky & Lovette 2008).
}
\usage{
plotRate(bt, pars)
}
\arguments{
  \item{bt}{ A numeric vector of branching times }
  \item{pars}{ a list of parameters, with named components lam0, k, mu0, and z }
}

\details{
	Plots speciation and extinction curves through time for models fitted using \link{fitSPVAR}, \link{fitEXVAR}, and \link{fitBOTHVAR}.
	If you want to plot the rate from SPVAR, note that this corresponds to assuming a very large z parameter; Try suppling z = 10000 or similar.
	Likewise, if you want to plot the fit from EXVAR, supply k = 0.
}
\value{
	Just plots the rates.
  
}
\references{ 



}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{
  \code{\link{fitSPVAR}}, \code{\link{fitEXVAR}}, \code{\link{fitBOTHVAR}}  
}
\keyword{ models }