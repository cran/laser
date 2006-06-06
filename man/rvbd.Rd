\name{Rate-Variable Birth Death}
\alias{rvbd}

\title{ Fit Rate-Variable Diversification Model With Extinction to Branching Times }
\description{
  Fits a two-rate variant of the birth death model to branching times derived from 
  phylogenetic data.  The model assumes that a clade has diversified under a net 
  diversification rate \emph{r1} until some time \emph{st}, at which point the net
  diversification rate shifts to \emph{r2}.  The extinction fraction \emph{a} is assumed
  constant and can take values from [0, 1).  The shift time \emph{st} giving the maximum log-likelihood
  is found by sequentially optimizing \emph{r1} and \emph{r2} across a large number of possible
  shift times.     
}
\usage{
rvbd(x, ai = c(0.1, 0.5, 0.95), ints = NULL, verbose = FALSE, file = "out_rvbd.txt")
}
\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{ai}{ a vector of initial \emph{a} values for optimization of the likelihood function }
  \item{ints}{ the number of intervals.  See details }
  \item{verbose}{ if 'verbose = TRUE', writes likelihoods and parameter
  estimates for all shift points considered to the indicated file. }
  \item{file}{ a filename for output if 'verbose = TRUE' }
}
\details{
  Non-linear optimization can be notoriously difficult, and the functions used in fitting 
  this model can be susceptible to the presence of multiple optima in the likelihood surface.  
  To guard against the possibility, the model is fitted using 3 initial values of the extinction
  fraction, \emph{a}.  You should check results obtained against those obtained using the 
  \code{\link{yule2rate}} model.  If the maximum log-likelihood under \code{\link{rvbd}} is less than 
  that obtained using \code{\link{yule2rate}}, your analysis has been fooled by a local (sub)optimum.  
  Repeat the analysis using
  additional starting values for \emph{a}, e.g., \samp{'ai = seq(0, .99. length.out = 15)'}. 
  
  'ints' is used in determining the number of shift points to consider.  If 'ints = NULL' (the 
  default), the model will consider only observed branching times as possible shift points.  See
  \code{\link{yule-n-rate}} for additional discussion of the 'ints' option.
  
  'verbose = TRUE' will write maximum log-likelihoods and parameter estimates for each shift time 
  under consideration to .  The
  file can then be loaded to examine the likelihood of a rate shift at different points in time.
    
}
\value{
  a list with the following components:
  \item{LH }{ the log-likelihood at the maximum}
  \item{aic }{ the Akaike Information Criterion}
  \item{r1 }{ the initial net diversification rate}
  \item{r2 }{ the final net diversification rate}
  \item{a }{ the extinction fraction}
  \item{st }{ the shift time giving the maximum log-likelihood}
}
\references{ 
Barraclough, T. G., and A. P. Vogler.  2002.  Recent diversification 
rates in North American tiger beetles estimated from a dated mtDNA 
phylogenetic tree. \emph{Mol. Biol. Evol.} 19:1706-1716.

Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed 
evolutionary process. \emph{Philos. Trans. R. Soc. Lond. B} 344:305-311.

Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 
}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ 
  Note that the net diversification rate \emph{r} is equal to the speciation rate \emph{S} minus 
  the extinction rate \emph{E}.  The extinction fraction \emph{a} = \emph{E/S}.  The 
  model assumes a constant extinction fraction \emph{a}, but this is not the same as
  a constant extinction rate \emph{E}.
  
  The model contains four free parameters: the initial and final net diversification rates,
  the extinction fraction, and the shift time.  

  Note that shift times, like branching times, are given in divergence units before present.  Thus, if 
  you have scaled a set of branching times to a basal divergence of 30 million years before present, you 
  would interpret 'st1 = 19.5' as an inferred shift point 19.5 million years before present.
  
  Large numbers of intervals or large sets of initial \emph{a} values can result in high computational time
  
}
\seealso{ \code{\link{yule-n-rate}}, \code{\link{fitdAICrc}} }
\examples{ 
  data(plethodon)
  result <- rvbd(plethodon, ints = 100, verbose = TRUE, file = 'temp.txt')
  
  ### now to generate a plot of log-likelihoods against shift times:
  temp <- read.table(file = 'temp.txt', header = TRUE)
  
  ### Rescaling shift times to reflect time since basal divergence
  stvec <- plethodon[1] - temp$st
  plot(temp$LH~stvec, xlab = "Time From Basal Divergence", 
        ylab = "Log-likelihood")
  
  unlink('temp.txt')  #clean-up temp file.      
  
}
\keyword{ models }
