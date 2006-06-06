\name{truncateTree}
\alias{truncateTree}

\title{ Truncate a Set of Branching Times }
\description{
  This function is used to truncate a set of branching times derived from a phylogenetic 
  tree.  Incomplete taxon sampling causes a spurious decline in the rate of lineage accumulation
  over time, and this effect becomes more severe towards the present.  Likewise, in the absence of 
  dense phylogeographic sampling, it may be desirable to omit the final one or more nodes from the tree.
  \samp{truncateTree} permits the user to omit \emph{k} final nodes, or the final \emph{t} time units.   
}
\usage{
truncateTree(x, omit.time = NULL, omit.nodes = NULL, batch = FALSE)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{omit.time}{ remove the final 'omit.time' time units before present from tree }
  \item{omit.nodes}{ remove the final 'omit.nodes' nodes before present }
  \item{batch}{ if TRUE, processes batch of branching times.  FALSE assumes
  as single set of branching times }
}
\details{
  Consider a set of branching times \samp{x = (100, 80, 50, 40, 30, 20, 10, 5)}.  
  If you wanted to analyze only the first half of the tree, perhaps due to concern about 
  incomplete sampling, calling \samp{truncateTree(x, omit.time = 50)} would return a vector
  of branching times \samp{x1 = (50, 30, 0)}.  Likewise, if you wished to omit the final 2 
  branching times, \samp{truncateTree(x, omit.nodes = 2)} would return 
  \samp{x2 = (90, 70, 40, 30, 20)}.  
  
  'batch' implies that you are entering a matrix or dataframe of branching times 
  for multiple trees, as returned by \code{\link{getBtimes.batch}}.  The 'omit.time'
  option is not available for batch processing.
}
\value{
  a numeric vector of branching times, or if \samp{batch = TRUE}, a matrix
  of branching times.
}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ 
  In the absence of dense phylogeographic sampling, it may be desirable to omit the final
  few nodes.
}
\examples{
  data(plethodon)
  pleth2 <- truncateTree(plethodon, omit.nodes = 2)
  #omits final 2 branching times
  
  plotLtt(pleth2)
}
\keyword{ manip }
