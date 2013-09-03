\name{getBtimes.batch}
\alias{getBtimes.batch}

\title{ Branching Times for Batch of Phylogenies }
\description{
  Processes files containing multiple sets of phylogenetic trees in 
  parenthetic ('Newick') format and returns a matrix of branching times to generate
  distributions of diversification rate test statistics.
}
\usage{
getBtimes.batch(fname, basal = NULL)
}

\arguments{
  \item{fname}{ filename where the parenthetic trees are stored }
 \item{basal}{ scales all trees to same basal divergence time.  See 'details'. }
}
\details{

  'basal' will scale all of your trees to the same basal divergence time.  This situation could arise
  in practice if you are interested in the posterior distributions of diversification rate parameters
  estimated under one or more models in this package.  For example, you might have the output file from a 
  run of MrBayes (the \emph{.t} file), with trees generated under a clock constraint.  Suppose you wished
  to examine the posterior distribution of speciation and/or extinction rates under a constant rate 
  birth death model.  Since all of the trees in the posterior distribution should be calibrated to the 
  same basal divergence, you may specify \samp{'basal = value'}, where value is the inferred time of the basal
  divergence.    
    
}
\value{
  
  A matrix of branching times, where rows are different trees or datasets, and columns are branching times.
  Thus, if you have N trees and K taxa, you will have a matrix of N rows and K-1 columns, since the number
  of branching times is one less than the number of taxa in a phylogenetic tree.  
  
  Suppose you have a file 'trees.tre', in newick format.  \samp{Res <- getBtimes.batch(fname = 'trees.tre')}
  returns the matrix of branching times.  You can access the j'th tree as \samp{Res[j, ]}.  Thus,
  \samp{plotLtt(Res[5,])} would generate a log-lineages through time plot for the 5th tree in the file.
}

\author{ Dan Rabosky \email{drabosky@umich.edu} }
\note{ 
  
  To analyze a single tree, see \code{\link{getBtimes}}
  
}
\seealso{
  \code{\link{getBtimes}}
}
\examples{

  data("rtrees50")
  write.table(rtrees50, file = 'temp.txt', quote=FALSE, 
                row.names = FALSE, col.names = FALSE)
  
  #creates a temporary file with trees in Newick format, identical to
  # output from PHYLOGEN & other software
  
  btimes <- getBtimes.batch(fname = "temp.txt")  
  
  # now btimes is a matrix of branching times.  Rows are different trees; 
  # columns are branching times.
  # To verify that this has correctly read the tree, we can plot the log-
  # lineages through time for the first tree:
  
  plotLtt(btimes[1,])
  
  # And we can compute the gamma statistic for this set of branching times:
  gamStat(btimes[1,])
  
  # or if you wanted to compute the gamma statistic for each tree:
  result <- as.numeric(apply(btimes, MARGIN=1, gamStat, return.list=FALSE))
  hist(result)  #plot histogram of gamma stat values for trees
  
  unlink("temp.txt") #clean up; delete temp file.
  
}
\keyword{ IO }
\keyword{ manip }
