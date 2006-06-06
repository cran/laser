\name{pureBirth}
\alias{pureBirth}

\title{ pureBirth }
\description{
  Fits pure birth (Yule) model to set of branching times
}
\usage{
pureBirth(x)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
}

\value{
  A list containing the following elements:
  \item{LH }{ the log-likelihood at the maximum}
  \item{aic }{ the Akaike Information Criterion}
  \item{r1 }{ the speciation rate giving the maximum log-likelihood}

}
\references{ 

Nee, S.  2001.  Inferring speciation rates from phylogenies. \emph{Evolution} 55:661-668.

Yule, G. U. 1924. A mathematical theory of evolution based on the conclusions 
of Dr.  J. C. Willis. \emph{Phil. Trans. R. Soc. Lond. B} 213:21-87.

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{yule-n-rate}}, \code{\link{bd}}, \code{\link{yuleWindow}}, \code{\link{rvbd}} }
\examples{
  data("plethodon")
  ### loads branching times for plethodontid salamander dataset
  pureBirth(plethodon)
}
\keyword{ models }