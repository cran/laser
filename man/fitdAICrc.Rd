\name{fitdAICrc}
\alias{fitdAICrc}

\title{ Test for Rate Variation Using delta-AICrc Test Statistic }
\description{
  Fits a specified set of rate-variable and rate-constant variants of the birth-death
  model to branching times from phylogenetic data.  The test statistic dAICrc is the difference
  in AIC scores between the best rate-constant and rate-variable models.
}
\usage{
fitdAICrc(x, modelset = c("pureBirth", "bd", "DDL", "DDX", "yule2rate", 
                "rvbd"), ints = NULL)
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
  models.  Available options are \code{\link{DDX}}, \code{\link{DDL}}, \code{\link{yule2rate}}, \code{\link{rvbd}}, and \code{\link{yule3rate}}. 
  See full descriptions of each of these models this document.
  
  'ints' is used in determining the number of shift points to consider.  If 'ints = NULL' (the 
  default), the model will consider only observed branching times as possible shift points.  See
  \code{\link{yule-n-rate}} for additional discussion of the 'ints' option.   
  
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
The dAICrc statistic described in Rabosky (2006) is one of the most powerful 
tests of temporal variation in diversification rates and may be the only statistical
approach that can infer temporal increases in diversification rates.  The approach involves fitting both 
rate-constant and rate-variable models to a set of branching times.  The statistic dAICrc is the 
difference in AIC scores between the best-fit rate-constant and best-fit rate-variable models.  

However, several considerations make this approach somewhat computationally intensive.  First, use of the 
AIC (as well as other model selection criteria, such as the AICc and BIC) result in model overfitting.  Simply
selecting the model with the lowest AIC criterion results in unacceptably high Type I error rates.  Moreover,
Type I error rates show a dependency on both the size of the phylogeny under consideration, as well as the 
number and type of rate-variable models in the candidate set.  It is thus necessary to derive the critical values of the distribution of dAICrc through simulation.      

The approach consists of 3 steps:

1. Analyze the sample phylogeny of N taxa using an appropriate candidate set of rate-variable models
using the function \code{fitdAICrc}.  Both the \code{pureBirth} and \code{bd} rate-constant models
must be included in the analysis.

2. Using the program \samp{Phyl-O-Gen}, the function \code{\link{yuleSim}}, or other phylogenetic
simulation software, generate a set of phylogenetic trees or branching times under the pure birth (Yule)
model.  Simulated trees should have the same number of tips as the test phylogeny, and should have the same
level of sampling completeness as the test phylogeny.  

3. Generate a matrix of branching times using \code{\link{getBtimes.batch}} if necessary.  Then run
\code{\link{fitdAICrc.batch}} on the set of simulated branching times, using the same set of candidate models
and 'ints' used to analyze the test phylogeny.  \code{\link{fitdAICrc.batch}} will generate the null
distribution of the dAICrc test statistic under the null hypothesis of rate-constancy and will determine 
the p-value of the dAICrc index computed for the test phylogeny.

Use of the 'ints' option is recommended.  Although the log-likelihood under the discrete-shift rate-variable
models will increase with greater numbers of ints, use of large numbers of ints (>200) typically result in very small improvements.  You 
can explore the effects of different numbers of 'ints' using \code{\link{rvbd}} or the \code{\link{yule-n-rate}} 
set of functions.  It is recommended that you explore the data in this way before undertaking this analysis; this
will also give you a sense of the computational time required to fit a given model with different numbers of 'ints'.

In practice, because the significance of dAICrc is assessed through simulation, the actual number of 'ints' is
not critical.  However, it is of critical importance that you use the same number of 'ints' when analyzing the 
test phylogeny (using \code{\link{fitdAICrc}}) and the simulated phylogenies (using \code{\link{fitdAICrc.batch}}).  

}
\seealso{ \code{\link{fitdAICrc.batch}} }
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

#we would still need to generate the null distribution of dAICrc
#through simulation using fitdAICrc.batch

}
\keyword{ models }
\keyword{ htest }
