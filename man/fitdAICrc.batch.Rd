\name{fitdAICrc.batch}
\alias{fitdAICrc.batch}
\alias{dAICrc.batch}
\alias{summarydAICrc}

\title{ Generate Null Distribution of dAICrc}
\description{
  Generates null distribution of dAICrc test statistic and determines significance
  of dAICrc indices calculated for 'real' phylogenies.  
}
\usage{
fitdAICrc.batch(x, modelset = c("pureBirth", "bd", "DDL", 
                "DDX", "yule2rate"), verbose = TRUE, 
                file = "out_daic.batch.txt", ints = NULL, 
                alternative = NULL, stat = NULL)
}

\arguments{
  \item{x}{ a numeric matrix or dataframe of branching times }
  \item{modelset}{ the set of rate-constant and rate-variable candidate models to be fitted }
  \item{verbose}{ write output to file (recommended) }
  \item{file}{ a filename for output if 'verbose = TRUE' }
  \item{ints}{ the number of shift points to consider }
  \item{alternative}{ specifies one or two-tailed Ho test }
  \item{stat}{ the dAICrc statistic for the 'real' phylogeny }
}
\details{
  For details on input format, see \code{\link{getBtimes.batch}}.  The \samp{modelset} should 
  be identical to that used to compute the dAICrc for the test phylogeny.  
  
  'verbose' generates an output file, specified by 'file', and contains parameter and likelihood summaries for
  each set of branching times.  It is highly recommended that you leave \samp{verbose = TRUE}, unless you have
  a reason not to.
  
  'alternative' is only applicable if you are entering a calculated dAICrc value for a test phylogeny (see
  'stat' below.  The default ('alternative = NULL') implements a two-tailed test of the null hypothesis that 
  rates have not changed over time.  You can also specify \samp{alternative = 'decrease'} or \samp{alternative = 
  'increase'} if you want a one-tailed test of the null hypothesis that rates have not decreased (or increased)
  over time.   
  
  'stat' is the calculated dAICrc statistic obtained for a test phylogeny for which you would like to
  obtain a p-value (using \code{\link{fitdAICrc}}).
  
}
\value{
  A dataframe with N rows, where N is the number of sets of branching times analyzed.  Columns of the data
  frame contain the following elements (and are accessed using the names given below):
  \item{dAICrc }{the calculated dAICrc statistic for dataset i}
  \item{rcbest }{the AIC score for the best rate-constant model}
  \item{rvbest }{the AIC score for the best rate-variable model}
  \item{bestmodel }{the best-fit model as determined by the AIC}
  \item{sp }{general diversification pattern specified by the best-fit model.  'NC' is \emph{no change
  in rate}, 'RD' is \emph{rate decrease}, and 'RI' is \emph{rate-increase}}. 
  \item{LpureBirth }{the maximum log-likelihood under the pure birth model}
  \item{Lbd }{the maximum log-likelihood under the rate-constant birth-death model}
  \item{L<model>}{Other model log-likelihoods abbreviated as above...}
}
\references{ 
Rabosky, D. L. 2006.  Likelihood methods for inferring temporal shifts in 
diversification rates.  \emph{Evolution} 60:1152-1164. 

}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }
\note{ 
Computational time required to generate the null distribution of dAICrc can be high 
if the number of phylogenies is large (>5000) and if \code{rvbd} and/or \code{yule3rate} models
are included in the candidate set. 

Ensure that your 'ints' settings are identical to those used to generate the dAICrc 
statistic for the test phylogeny (\code{\link{fitdAICrc}}).

\code{fitdAICrc.batch} will print details of run progress to screen, but you
must ensure that you turn off the 'buffered output' option under the 'Misc' 
pull-down menu in the R GUI.
}
\seealso{ \code{\link{fitdAICrc}} }
\examples{

data(rtrees50)
write.table(rtrees50, file = 'temp.txt', quote=FALSE, 
                row.names = FALSE, col.names = FALSE)

bt <- getBtimes.batch(file = 'temp.txt', format = 'newick')
# using 'string' option to get branching times.  bt is now a dataframe
# of branching times, where rows are trees.

result <- fitdAICrc.batch(bt, modelset = c("pureBirth", "bd", "DDX", 
                          "DDL", "yule2rate"), ints = 20)
# fits 2 rate-constant and 3 rate variable models to each of the
# 50 trees in the example datafile 'rtrees50' and generates 
# the distribution of dAICrc from the sample.
}
\keyword{ models }
\keyword{ htest }
