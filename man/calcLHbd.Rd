\name{calcLHbd}
\alias{calcLHbd}

\title{ Calculate Likelihood of Branching Times Under Birth-Death Model }
\description{
  Calculates the likelihood of branching times given a vector of branching times
  and parameters \emph{r} (net diversification rate; speciation rate \emph{S} - 
  extinction rate \emph{E}) and \emph{a} (the extinction fraction, \emph{E}/\emph{S})
}
\usage{
calcLHbd(x, r, a)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{r}{ the net diverrsification rate, \emph{S} - \emph{E} }
  \item{a}{ the extinction fraction, \emph{E}/\emph{S} }
}
\details{
  A function that can be called to explore alternative parameterizations of the
  birth-death process, to bootstrap likelihood confidence regions, or to generate
  plots of likelihood surfaces (as in Nee et al. 1994).
}
\value{
  the likelihood of the branching times given \emph{a} and \emph{r}
}
\references{ 
Kendall, D. G.  1948.  On the generalized "birth-and-death" process. \emph{Ann. Math. Stat.} 19:1-15.

Nee, S., E. C. Holmes, R. M. May, and P. H. Harvey.  1994a.  Extinction rates can be estimated from molecular phylogenies. \emph{Philos. Trans. R. Soc. Lond. B} 344:77-82.

Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed evolutionary process. \emph{Philos. Trans. R. Soc. Lond. B} 344:305-311.

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{bd}, \code{rvbd}, \code{pureBirth} }
\examples{

#plot a likelihood surface for the plethodon dataset
data("plethodon")
rvector <- seq(.001, .04, length.out = 100)
avector <- seq(0, 0.99, length.out = 100)
#calculating likelihoods:
lmatrix <- matrix(0, length(rvector), length(avector))
for (i in 1:length(rvector)){
  for (j in 1:length(avector)){
    lmatrix[i,j] <- calcLHbd(plethodon, rvector[i], avector[j])
    }
   }
 
 #Now to plot surface:
 lmax <- max(lmatrix)      #maximum calculated LH
 filled.contour(rvector, avector, lmatrix, 
                levels = seq(lmax-5, lmax, length.out = 20), 
                col = heat.colors(20), xlab = "Net Diversification Rate", 
                ylab = "Extinction Fraction", key.title = title("Log-LH"))
 
 # plots the surface.  See ?filled.contour for help on this plotting function.
               
}    

\keyword{ models }
