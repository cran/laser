\name{laser-package}
\alias{laser-package}
\alias{laser}
\docType{package}
\title{
LASER
}
\description{
	Maximum likelihood-based methods for analyzing lineage diversification rates
}
\details{
\tabular{ll}{
Package: \tab laser\cr
Type: \tab Package\cr
Version: \tab 2.4-0\cr
Date: \tab 2013-05-05\cr
License: \tab GPL 2.0 \cr
LazyLoad: \tab yes\cr
}

}
\author{
Dan Rabosky
Klaus Schliep
Maintainer: Dan Rabosky <drabosky@umich.edu>
}


\keyword{ package }
%\seealso{}
\examples{
	data(warblers)
	fitSPVAR(warblers)
	fitEXVAR(warblers)
	gamStat(warblers)
	mccrTest(CladeSize=50, NumberMissing=10, NumberOfReps=50, ObservedGamma=-2.1)
}
