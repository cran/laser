\name{MagallonSanderson.Rd}
\alias{lambda.crown.ms01}
\alias{lambda.stem.ms01}
\alias{lambda.stem.ml}
\alias{lambda.stem.ci}

\title{ Estimate speciation rates for crown and stem groups }
\description{
  Esimate speciation rates under any relative extinction rate, given a crown or stem group age and a number of extant species
}    
\usage{
lambda.crown.ms01(n, tb, eps=0)
lambda.stem.ms01(n, tb, eps=0)
lambda.stem.ml(n, tb, eps=0)
lambda.stem.ci(tb, r, eps=0)
}


\arguments{
  \item{n}{ the number of extant species in a clade}
  \item{tb}{ the stem or crown group age}
  \item{eps}{ the relative extinction rate, or \code{mu / lambda} }
  \item{r}{ the net diversification rate, r (\code{lambda.stem.ci} only)}
}
\details{
\code{lambda.crown.ms01} estimates speciation rates assuming some value of \code{eps} and a crown group age after Magallon & Sanderson (2001), eqn A7.

\code{lambda.stem.ms01} estimates speciation rates assuming a known stem clade age.  Same as eqn A6 in Magallon & Sanderson 2001.

\code{lambda.stem.ml} finds maximum likelihood estimate of speciation rate given stem age and \code{eps}. This is given in Magallon & Sanderson (2001), eqns 1-2,  Raup (1985) eqn A18, and others.

\code{lambda.stem.ci} gives 95 percent confidence intervals on expected species diversity after \code{tb} time units given net diversification rate \code{r} and \code{eps}

}
\value{
a list with the following components (for speciation rate estimators)
\item{lambda}{the estimated speciation rate}
\item{r}{the net diversification rate}

or, for \code{lambda.stem.ci}
\item{upper}{upper bound of 95 percent confidence interval on expected species diversity}
\item{lower}{lower bound of 95 percent confidence interval on expected diversity}

}
\references{ 
Magallon, S., and M. J. Sanderson. 2001. Absolute diversification rates in angiosperm clades. Evolution 55:1762-1780.

Raup, D. M. 1985. Mathematical models of cladogenesis. Paleobiology 11:42-52.

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_2rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}}, \code{\link{fitNDR_1rate}}}
\examples{

}
\keyword{ models }
