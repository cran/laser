\name{fitdAICrc}
\alias{fitdAICrc}

\title{ Test for Rate Variation Using delta-AICrc Test Statistic }
\description{
  Fits a specified set of rate-variable and rate-constant variants of the birth-death
  model to branching times from phylogenetic data.  The test statistic dAICrc is the difference
  in AIC scores between the best rate-constant and rate-variable models.
}
\usage{
fitdAICrc(x, modelset = c("pureBirth", "bd", "DDL", "DDX", "yule2rate"), ints = NULL)
}
\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{modelset}{ the set of rate-constant and rate-variable candidate models to be fitted }
  \item{ints}{ the number of intervals.  See 'Details' }
}
\details{
  \samp{fitdAICrc} implements the dAICrc test statistic for temporal variation in diversification 
  rates as described in Rabosky (2006).  
  
  \code{modelset} is a list of the rate-constant and rate-variable models to consider.  You should include 
  both rate-constant models (\code{\link{pureBirth}} and \code{\link{bd}}), as well as one or more candidate rate-variable
  models.  Available options are \code{\link{DDX}}, \code{\link{DDL}}, \code{\link{yule2rate}}, and \code{\link{yule3rate}}. 
  See full descriptions of each of these models this document.
  
  'ints' is used in determining the number of shift points to consider.  If 'ints = NULL' (the 
  default), the model will consider only observed branching times as possible shift points.  See
  \code{\link{yule-n-rate}} for additional discussion of the 'ints' option.  ]
  
  Note that the rvbd function is no longer suppored ('rate variable birth death').
  
}
\value{
  a dataframe with the number of rows equal to the number of candidate models.  Columns include likelihoods,
  parameters, and AIC scores for each model.  The first column contains the model names.  If a parameter
  is not present in a particular model, it will have an entry of 'NA' in the column for that parameter.  
  Parameter names follow conventions for model descriptions in other parts of this document.  For example, 
  parameter \samp{r1} is the initial net diversification rate for all models (note that this will be the 
  only rate for the \samp{pureBirth} model).  
  
  The full set of columns if all available models are included in the candidate set will consist 
  of the following:
  
  \item{model}{ the model name for row \emph{i} in the dataframe}
  \item{params}{ the free parameters for model[i] }
  \item{np}{ the number of free parameters in mode[i]}
  \item{mtype}{ either 'RC' for rate-constant or 'RV' for rate-variable}
  \item{LH}{ the log-likelihood under model[i]}
  \item{r1, r2, r3}{net diversification rates, as applicable; r1 is always the initial rate, and 
  r3 is always the final rate}
  \item{a}{ the extinction fraction \emph{E/S} if applicable}
  \item{xp}{the x-parameter from the \code{\link{DDX}} model}
  \item{k}{the k-parameter from the \code{\link{DDL}} model}
  \item{st1, st2}{ shift-times, if applicable.  st1 is always the first shift point}
  \item{AIC}{the Akaike Information Criterion for model[i]}
  \item{dAIC}{delta-AIC; the difference in AIC scores between model[i] and the overall best-fit model}
}
\references{ 

Nee, S., R. M. May, and P. H. Harvey.  1994b.  The reconstructed 
evolutionary process. \emph{Philos. Trans. R. Soc. Lond. B} 344:305-311.

Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\note{ 


}
\seealso{ }
\examples{
data(agamids)
agbtimes <- getBtimes(string = agamids)
#agbtimes is now a vector of branching times from the agamid phylogeny

#here we fit 2 rate-constant and 3 rate-variable models 
# to the agamid data:
result <- fitdAICrc(agbtimes, modelset = c("pureBirth", "bd",
          "DDX", "DDL", "yule2rate"), ints = 100)

# this outputs summaries of parameters and likelihoods to screen; 
# object 'result' is a dataframe containing all parameter estimates, 
# likelihoods, and AIC scores



}
\keyword{ models }
\keyword{ htest }
