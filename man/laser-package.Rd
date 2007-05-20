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
Version: \tab 2.0\cr
Date: \tab 2008-05-19\cr
License: \tab GPL 2.0 ? \cr
LazyLoad: \tab yes\cr
}

}
\author{
	Dan Rabosky
Maintainer: Dan Rabosky <DLR32@cornell.edu>

}


\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
	data(warblers)
	fitSPVAR(warblers)
	fitEXVAR(warblers)
	gamStat(warblers)
	mccrTest(CladeSize=50, NumberMissing=10, NumberOfReps=50, ObservedGamma=-2.1)
}
