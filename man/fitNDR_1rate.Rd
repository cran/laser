\name{fitNDR_1rate.Rd}
\alias{fitNDR_1rate}

\title{ Estimate diversification rate from combined taxonomic/phylogenetic data  }
\description{
  Find maximum likelihood estimate of net diversification rate from combined taxonomic/phylogenetic data.
}    
\usage{
fitNDR_1rate(phy, eps=0, rbounds=c(0.0001, .5), combined=TRUE)

}


\arguments{
  \item{phy}{ a class 'phylo' phylogenetic tree with an additional component 'phenotype'; \code{phy$phenotype} specifies the number of tips per terminal taxon.  The phenotype data is added with \code{getTipdata}}
  \item{eps}{ the relative extinction rate, or \code{mu / lambda} }
  \item{rbounds}{ The upper and lower limits for the 1-dimensional optimization of the net diversification rate }
  \item{ combined }{ Should likelihoods be calculated using both combined phylogenetic and taxonomic data, or taxonomic data only? }
}
\details{
  \code{fitNDR_1rate} finds the maximum likelihood estimate of the net diversification rate using phylogenetic and taxonomic data. This method is best applied to phylogenetic trees for higher taxonomic levels where you have incomplete sampling but know the approximate species diversities of each terminal taxon in the tree (e.g., a phylogenetic tree of arthropod families).
  
  \code{eps} allows you to estimate the speciation rate under any assumed (constant) relative extinction rate, where the relative extinction rate is \code{mu / lambda }. 
  
  \code{rbounds} is an important argument, and the default may not work. Basically, you are trying to find the ML estimate of the net diversification rate; however, 
  the optimization algorithm requires that you specify min and max values for the search-space. Sometimes the function cannot be evaluated at these bounds or the optimization will fail for other reasons.
  Thus, you might wish to start with a vary small range (e.g., rbounds = c(0.001, 0.05)) and increase the range until you are confident of finding the optimum.
  Most important: if the maximum diversification rate appears to be on a boundary, you have a problem. For example, if you specify rbounds = c(0.001, 0.05), and find the 
  \code{r} value returned by the function to be 0.05 at the maximum, then this is almost certainly not valid. You will need to repeat your analysis, expanding the range as necessary, until you find a maximum that does
  not lie on a boundary. Also be wary of multiple optima, though I haven't personally encountered them.
  
  
  The \code{combined} argument asks whether you wish to use an estimator based on (i) combined taxonomic and phylogenetic data, or (ii) taxonomic data only. There are good reasons for trying both approaches: when \code{combined = TRUE}, both the taxonomic/species richness and phylogenetic backbone of a tree contribute heavily to the overall likelihood estimate.  This occurs in spite of the fact that most of the species diversity is actually in the tips. When \code{combined = TRUE}, the internal phylogenetic structure of the tree does not contribute to the likelihood; you are finding the maximum likelihood estimate of the speciation rate \code{lambda} from the species richness data plus stem clade ages alone. I recommend checking estimates under both \code{combined =TRUE} and \code{combined=FALSE}

}
\value{
  a dataframe with the following components:
  \item{LH }{ The log-likelihood at the maximum}
  \item{aic }{the Akaike Information Criterion}
  \item{r }{ the net diversification rate,\code{lambda - mu}, at the maximum}
  \item{lambda }{ the ML estimate of the speciation rate}
  \item{eps }{ the value of \code{eps} used}
  
}
\references{ 

Rabosky, D. L., S. C. Donnellan, A. L. Talaba, and I. J. Lovette. 2007. Exceptional among-lineage variation in diversification rates during the radiation of Australia's largest vertebrate clade. Proc. Roy. Soc. Lond. Ser. B 274:2915-2923.


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_2rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}} }
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

	# and you can see that the 'best' rate shift 
	# location is the MRCA of the genera Ctenotus and Lerista
}
\keyword{ models }
