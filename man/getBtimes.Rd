\name{getBtimes}
\alias{getBtimes}
\title{ Get Branching Times From Parenthetic-Format Tree ('Newick') }
\description{
  Reads a file containing an ultrametric tree in parenthetic format and 
  returns a numeric vector of branching times, sorted from earliest to most 
  recent.
}
\usage{
getBtimes(string = NULL, file = NULL)
}

\arguments{
  \item{string}{  the name of an object containing a 'newick' tree character string}
  \item{file}{  a file containing a single tree in 'newick' (parenthetic) format }
}
\details{
  This function requires the package \samp{ape}.  If \samp{ape} has not been 
  attached to the search path but exists in the R directory, it will automatically
  be loaded.  If \samp{ape} does not exist, you must install it by typing
  \samp{install.packages("ape")} at the R prompt.
  
  If tree is not ultrametric, an error message will be generated.  
  
}
\value{
  a numeric vector of branching times
}


\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{getBtimes.batch}} }
\examples{
data(agamids)
write.table(agamids, file = 'example.tre', quote=FALSE, 
                row.names = FALSE, col.names = FALSE)
 # the preceding lines generate a file 'example.tre' that can be 
 # read by 'getBtimes'
 
 agbtimes <- getBtimes(file = 'example.tre')
 #or alternatively:
 agbtimes <- getBtimes(string = agamids)
 #agbtimes is now a numeric vector of branching times
 
 plotLtt(agbtimes) #plot log-lineages through time for agamid data
 gamStat(agbtimes) #calculate gamma statistic for agamid data 
 unlink('example.tre') #clean-up. 

}
\keyword{ IO }
\keyword{ manip }
