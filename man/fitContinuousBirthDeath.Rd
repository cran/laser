\name{fitSPVAR}
\alias{fitSPVAR}
\alias{fitEXVAR}
\alias{fitBOTHVAR}

\title{ Fit model with continuous-time varying speciation/extinction rates to phylogeny }
\description{
  Functions that fit time-varying speciation and extinction models to branching times
  derived from phylogenetic data. \samp{fitSPVAR} fits a model with an exponentially
  declining speciation rate through time and constant extinction. \samp{fitEXVAR} fits a model 
  with exponentially increasing extinction and constant speciation. \samp{fitBOTHVAR} fits
  a model where both speciation and extinction rates can vary through time. 
  
}    
\usage{
fitSPVAR(bt, init=c(2, .2, .1))
fitEXVAR(bt, init=c(.3, .01, 1))
fitBOTHVAR(bt, init=c(.3, .5, .1, .5))
}

\arguments{
  \item{bt}{ a numeric vector of branching times }
  \item{init}{ a vector of starting parameter values for the optimization algorithm. Note that there are 3 parameters in SPVAR and EXVAR models, but 4 parameters in BOTHVAR}
}

\details{
  These functions are described in Rabosky & Lovette (2008, \emph{Evolution}). 
  There is a high likelihood that the optimization will fail across a large span of parameter 
  space, so you should not expect the default parameters to work. I recommend varying the first parameter (e.g., init[1])
  before changing anything else.  
  
  The three models return the log-likelihood, the AIC, and the parameter estimates. Parameters correspond exactly to those 
  described in Rabosky and Lovette (2008), equations 7-11. The speciation rate is specified by parameters lam0 and k, and extinction
  through time is described by mu0 and z.  lam0 and mu0 are the initial speciation and final extinction rates, respectively. 
  k and z control the rate of decrease / increase in speciation and extinction, respectively. 
  
  Because the underlying mathematical model of speciation and extinction requires that the extinction rate not exceed 
  the speciation rate, the models have been reparameterized as follows: for the SPVAR model, optimization is performed on 
  the parameter set c(X, k, mu0), where k and mu0 are defined as above and X corresponds the net diversification rate. This provides for a
  fixed minimum bound for the optimization interval (because X must always be greater than 0).
  This should be irrelevant to the useage of the function. However, if you are having a difficult time getting optimization 
  to work without an error, bear in mind that the initial parameters you are suppling correspond to:
  SPVAR: X, k, mu0; 
  EXVAR: X, mu0, z;
  and BOTHVAR, X, k, mu0, z.
  A good solution might be to start with initial parameter values under the pure birth or constant rate birth death model. You will 
  need to estimate these in light of equations 7-11 in Rabosky & Lovette (2008).
  
  
  Optimization uses the L-BFGS-B algorithm as implemented in the function \code{optimize}. Numerical integration uses \code{integrate}.  
  
  You can use the equations from Rabosky & Lovette (2008) with your parameter estimates to generate a plot of the speciation and extinction
  through time curves, or you can use the function \link{plotRate}
  
}
\value{

  \item{model }{  The name of the model}
  \item{LH }{ Maximum log-likelihood of the fitted model}
  \item{aic }{ AIC }
  \item{ lam0 }{ The initial speciation rate }
  \item{ k }{ Parameter of the exponential change in speciation rate }
  \item{ mu0 }{ the final extinction rate }
  \item{ z }{ Parameter of the exponential change in extinction rate }
  
}
\references{ 

Rabosky, D. L. and I. J. Lovette. 2008. Explosive evolutionary radiations: decreasing speciation or increasing extinction through time? Evolution, doi:10.1111/j.1558-5646.2008.00409.x

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\examples{
	data(warblers)
	fitSPVAR(warblers)
	fitEXVAR(warblers)
	fitBOTHVAR(warblers)
}
\keyword{ models }
