\name{fitNDR_2rate.Rd}
\alias{fitNDR_2rate}

\title{ Fit 2-rate diversification model to combined phylogenetic/taxonomic data  }
\description{
  Finds likelihoods and parameter estimates for combined phylogenetic/taxonomic data under a two-rate diversification model. 
}    
\usage{
fitNDR_2rate(phy, eps=0, combined=TRUE, rate.decrease=FALSE, rbounds=c(.001, .5))

}


\arguments{
  \item{phy}{ a class 'phylo' phylogenetic tree with an additional component 'phenotype'; \code{phy$phenotype} specifies the number of tips per terminal taxon.  The phenotype data is added with \code{getTipData}}
  \item{eps}{ the relative extinction rate, or \code{mu / lambda} }
  \item{ combined }{ Should likelihoods be calculated using both combined phylogenetic and taxonomic data, or taxonomic data only? }
  \item{rate.decrease}{Constrain the analysis such that the highest diversification rate occurs in the tree bipartition containing the root node?}
  \item{rbounds}{ The upper and lower limits for the 1-dimensional optimization of the net diversification rate. See details. }

}
\details{
  \code{fitNDR_2rate} fits a 2-speciation rate model to combined phylogenetic/taxonomic data. Suppose you have a higher-level phylogeny for some group of organisms (e.g., beetle families), where you also know the approximate species diversities for each terminal. \code{fitNDR_2rate} assumes that at some point in the tree, an ancestral speciation rate \code{lambda1} shifts to a new speciation rate \code{lambda2}. 
  
  The model is fitted by iteratively splitting the tree at each node and fitting a birth-death model to each of the resulting bipartitions. Thus, for each node, you obtain (i) the likelihood of a rate shift at that position in the tree, and (ii) the estimated speciation rates for each bipartition. The function returns a dataframe giving the likelihoods of rate shift at each node as well as the parameter estimates. 
  
    \code{eps} allows you to estimate the speciation rate under any assumed (constant) relative extinction rate, where the relative extinction rate is \code{mu / lambda }. 
  
   \code{rbounds} is an important argument, and the default may not work. Basically, you are trying to find the ML estimate of the net diversification rate; however, 
  the optimization algorithm requires that you specify min and max values for the search-space. Sometimes the function cannot be evaluated at these bounds or the optimization will fail for other reasons.
  Thus, you might wish to start with a vary small range (e.g., rbounds = c(0.001, 0.05)) and increase the range until you are confident of finding the optimum.
  Most important: if the maximum diversification rate appears to be on a boundary, you have a problem. For example, if you specify rbounds = c(0.001, 0.05), and find the 
  \code{r} value returned by the function to be 0.05 at the maximum, then this is almost certainly not valid. You will need to repeat your analysis, expanding the range as necessary, until you find a maximum that does
  not lie on a boundary. Also be wary of multiple optima, though I haven't personally encountered them.
  
 
 Specifying \code{combined = TRUE} estimates \code{lambda} from combined taxonomic and phylogenetic data. Thus, the likelihood is a function of the internal branching structure of the tree and the species richness/taxonomic data. If \code{combined=FALSE}, likelihoods are estimated from the species richness/taxonomic data only, and the internal branching structure of the tree does not contribute.
  
  Option \code{rate.decrease} fits a model where the highest speciation rate must occur in the tree bipartition containing the root node.

}
\value{
  a dataframe with the following components:
  \item{node}{ The node defining the tree bipartition }
  \item{LH }{ The log-likelihood at the maximum}
  \item{aic }{the Akaike Information Criterion}
  \item{r.1}{ the net diversification rate,\code{lambda - mu}, for the partition containing the root node}
  \item{lambda.1}{the speciation rate for the for the partition containing the root node}
  \item{LH.1}{the log-likelihood for partition 1 at the maximum}
   \item{r.2}{ the net diversification rate,\code{lambda - mu}, for the partition NOT containing the root node, e.g., partition 2}
  \item{lambda.2}{the speciation rate for the for the partition 2} 
  \item{LH.2}{the log-likelihood for partition 2 at the maximum}
  \item{eps }{ the value of \code{eps} used}
  
}
\references{ 

Rabosky, D. L., S. C. Donnellan, A. L. Talaba, and I. J. Lovette. 2007. Exceptional among-lineage variation in diversification rates during the radiation of Australia's largest vertebrate clade. Proc. Roy. Soc. Lond. Ser. B 274:2915-2923.


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_1rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}}  }
\examples{
	data(skinktree);
	data(skinkdiversity);
	skinktree <- getTipdata(skinkdiversity, skinktree);
	
	#first we fit the one rate model
	fitNDR_1rate(skinktree, eps=0);

	#here we fit the 2 rate model
	res <- fitNDR_2rate(skinktree, eps=0);

	#extracting the node most likely to have undergone rate shift:
	subset(res, res$LH==max(res$LH));

	#this function plots the node numbers on the tree:
	plotNodeNumbers.phylo(skinktree);


}
\keyword{ models }