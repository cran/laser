\name{yule-n-rate}
\alias{yule-n-rate}
\alias{yule2rate}
\alias{yule3rate}
\alias{yule4rate}
\alias{yule5rate}


\title{ yule-n-rate }
\description{
  Fits multi-rate variants of the pure birth (Yule) model to branching times derived
  from phylogenetic data.  For example, the \emph{yule2rate} model assumes 
  that the clade has diversified under speciation rate \emph{r1} until some time \emph{st}, 
  at which point the speciation rate shifts to a new rate \emph{r2}.  The shift point(s) are 
  found by optimizing parameters and computing likelihoods for a set of possible shift times and selecting
  the parameter combinations giving the maximum log-likelihood.  
}
\usage{
yule2rate(x, verbose = FALSE, ints = NULL, file = "out_yule2rate.txt")

yule3rate(x, ints = NULL, verbose = FALSE, file = "out_yule3rate.txt") 

yule4rate(x, ints = NULL)

yule5rate(x, ints = NULL)
}

\arguments{
  \item{x}{ a numeric vector of branching times }
  \item{verbose}{ if 'verbose = TRUE', writes likelihoods and parameter
  estimates for all shift points considered to the specified file.  Default is 
  FALSE.  Only available for \emph{yule2rate} and \emph{yule3rate} models}
  \item{ints}{ the number of intervals.  See details }
  \item{file}{ a filename for output if 'verbose = TRUE' }
}
\details{
  
  'verbose' - for \emph{yule2rate} and \emph{yule3rate} models, maximum log-likelihoods
  and parameter estimates for each shift time under consideration will be output to file.  The
  file can then be loaded to examine the likelihood of a rate shift at different points in time.
  
  'ints' is used in determining the number of shift points to consider.  If 'ints = NULL' (the 
  default), the model will consider only observed branching times as possible shift points.  Suppose we have
  a small dataset with the following branching times: (100, 80, 50, 40, 30, 20, 10, 5, 2).  Under the \emph{yule2rate} model, we assume 
  that the clade has diversified under some rate \emph{r1} until some time \emph{ts}, at which point the rate
  simultaneously shifts in all lineages to a new rate \emph{r2}.  In this example, if 'ints = NULL', we would
  use the set of observed branching times only, but omitting the first and final branching times (thus, we 
  would be considering only \emph{st} = (80, 50, 40, 30, 20, 10, 5) as possible shift points.  
  If 'ints = 100', we would consider 100 evenly spaced shift points on the interval between the 2nd branching time
  and the N-1 branching time (e.g., on (80, 5)).  'ints' works well for \emph{yule2rate} and \emph{yule3rate}
  models, but can result in high computational times for \emph{yule4rate} and \emph{yule5rate} models.
}

\value{
  a data frame containing the following elements:
  
  \item{LH }{ the maximum log-likelihood}
  \item{AIC }{ the Akaike Information Criterion} 
  \item{r1}{ the first (earliest) speciation rate giving the maximum log-likelihood}
  \item{r2}{ the ML estimate of the second speciation rate (in the case of \emph{yule2rate}, this 
  will be the final speciation rate)}
  \item{st1}{ the earliest shift point (for \emph{yule2rate}, you will only have a single shift point)}
  \item{...}{ speciation rates and shift points for models other than \emph{yule2rate} are abbreviated by
  \emph{r3}, \emph{st2}, etc, as described above}
}
\references{ 
Barraclough, T. G., and A. P. Vogler.  2002.  Recent diversification rates in 
North American tiger beetles estimated from a dated mtDNA phylogenetic tree. 
\emph{Mol. Biol. Evol.} 19:1706-1716.

Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 
 }
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ The total number of parameters for each model is equal to the number of speciation rates
and shift points subject to optimization.  Thus, the \emph{yule3rate} model has 3 speciation rates and 2
shift times, for a total of 5 parameters.  Strictly speaking, it may be inappropriate to treat the 
shift time \emph{st} as a free parameter, as it can only take on a limited set of values.  However, 
in practice, it appears to work well; in many cases, using observed shift times 
can give higher likelihoods than when 'ints' are specified.  There seems to be little improvement in the
log-likelihood with 'ints' greater than 1000, at least for phylogenies with fewer than 100 tips.

Note that shift times, like branching times, are given in divergence units before present.  Thus, if 
you have scaled a set of branching times to a basal divergence of 30 million years before present, you 
would interpret 'st1 = 19.5' as an inferred shift point 19.5 million years before present.
}
\seealso{ \code{ bd }, \code{ fitdAICrc }, \code{ yuleWindow }, \code{ pureBirth } }
\examples{

data(plethodon)

### fitting a 3-rate Yule model to the Plethodon data:

result <- yule3rate(plethodon)

### gives data frame with maximum log-likelihood and parameter estimates 
### at the max.
### In this case, we would access individual parameters as
### result$LH (the max), result$st1 (first shift time), result$st2
### (the second shift time), and result$r1, result$r2, and result$r3
### for the speciation rates.

### Here we will use 'yule2rate' to output maximum log-likelihoods 
### for each shift point considered, then load the file and plot 
### log-likelihoods of a rate shift against 'time from basal divergence'
### to graphically explore the tempo of diversification

# result <- yule2rate(plethodon, ints = NULL, 
#          verbose = TRUE, file = 'out.txt')
# LHtable <- read.table(file = 'out.txt', header = TRUE)

### 'header = TRUE' ensures that variable names are correctly read

### rescaling shift times so that they reflect 'time from basal divergence':

# LHtable$st1 <- plethodon[1] - LHtable$st1
# plot(LHtable$LH~LHtable$st1, xlab = 'Time From Basal Divergence', 
#    ylab = 'Log-likelihood')


}
\keyword{ models }
