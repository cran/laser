\name{fitdAICrc.batch}
\alias{fitdAICrc.batch}


\title{ Generate Null Distribution of dAICrc}
\description{
  Generates null distribution of dAICrc test statistic and determines significance
  of dAICrc indices calculated for 'real' phylogenies.  
}
\usage{
fitdAICrc.batch(x, observed=NULL, ...)
}

\arguments{
  \item{x}{ a numeric matrix or dataframe of branching times }
  \item{observed}{ the dAICrc statistic for the 'real' phylogeny }
  \item{...}{ other parameters, e.g., ints for yule2rate }
}
\details{
  For details on input format, see \code{\link{getBtimes.batch}}.  
  
  This fits the following models: pure birth, birth death, density-dependent (exponential),
  density-dependent (logistic), yule-2-rate
  
  
  'observed' is the calculated dAICrc statistic obtained for a test phylogeny for which you would like to
  obtain a p-value (using \code{\link{fitdAICrc}}).
  
}
\value{
  A dataframe with N rows, where N is the number of sets of branching times analyzed.  Columns of the data frame
  are the AIC scores for each of the models. The final column, with name dAIC, is the deltaAICrc test statistic
  for that particular tree.

}
\references{ 
Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 

}
\author{ Dan Rabosky \email{drabosky@umich.edu} }
\note{ 
	Make sure that you are using the exact same set of models for the real phylogeny and for the simulated 
	phylogenies. If not, your results are invalid.

}

% \examples{}
\keyword{ models }
\keyword{ htest }
