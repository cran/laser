\name{getBtimes.batch}
\alias{getBtimes.batch}

\title{ Branching Times for Batch of Phylogenies }
\description{
  Processes files containing multiple sets of phylogenetic trees in 
  parenthetic ('Newick') format and returns a matrix of branching times to generate
  distributions of diversification rate test statistics.
}
\usage{
getBtimes.batch(file = NULL, format = "newick", basal = NULL)
}

\arguments{
  \item{file}{ filename where the parenthetic trees or set of branching times
  are stored }
  \item{format}{ format of the data you would like to process.  'newick' (default)
  assumes you have a file of parenthetic-format trees.  format 'branchingtimes' implies that
  data are a numerical list of branching times.  See 'details' belows.}
  \item{basal}{ scales all trees to same basal divergence time.  See 'details'. }
}
\details{
  Input formats: see example files (50RandomTrees.tre) for an example of 'newick' format.  
  This is the output format used by the phylogenetic simulation package 
  \emph{Phylogen}; you could also easily modify the \emph{.t} file of sampled trees from a run
  of \emph{MrBayes} or other MCMC sampler to satisfy the format requirements.
  
  It is recommended that you use the 'file' option, rather than 'string', for reading in large 
  numbers of trees.  See \code{\link{fitdAICrc.batch}} for an example of input using 'string'.
  
  See the file '50btimesyule.txt' for an example of the 'branchingtimes' format.  If you have
  'branchingtimes' formatted data, the first number in your file must be the number of taxa (NOT the 
  number of branching times); the second number must be the number of sets of branching times.  
  
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
  
  Suppose you have a file 'trees.tre', in newick format.  \samp{Res <- getBtimes.batch(file = 'trees.tre')}
  returns the matrix of branching times.  You can access the j'th tree as \samp{Res[j, ]}.  Thus,
  \samp{plotLtt(Res[5,])} would generate a log-lineages through time plot for the 5th tree in the file.
}

\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ 
  You \emph{MUST} have the package 'ape' installed to use this function, if format = 'Newick'.  If you 
  are connected to the internet, you can obtain 'ape' by typing \samp{'install.packages("ape")} at the R prompt.  
  If 'ape' is installed on your computer, \samp{getBtimes.batch} will automatically attach the package
  to the search path.
  
  To analyze a single tree, see \code{\link{getBtimes}}
  
}
\seealso{
  \code{\link{gamStat.batch}}, \code{\link{fitdAICrc.batch}}, \code{\link{getBtimes}}
}
\examples{

  # if format = 'newick', requires library(ape)
  data("rtrees50")
  write.table(rtrees50, file = 'temp.txt', quote=FALSE, 
                row.names = FALSE, col.names = FALSE)
  
  #creates a temporary file with trees in Newick format, identical to
  # output from PHYLOGEN & other software
  
  btimes <- getBtimes.batch(file = "temp.txt")  
  
  # now btimes is a matrix of branching times.  Rows are different trees; 
  # columns are branching times.
  # To verify that this has correctly read the tree, we can plot the log-
  # lineages through time for the first tree:
  
  plotLtt(btimes[1,])
  
  # And we can compute the gamma statistic for this set of branching times:
  gamStat(btimes[1,])
  
  # or if you wanted to compute the gamma statistic for each tree:
  result <- gamStat.batch(btimes)
  hist(result$gamstat)  #plot histogram of gamma stat values for trees
  
  unlink("temp.txt") #clean up; delete temp file.
  
}
\keyword{ IO }
\keyword{ manip }
