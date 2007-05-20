\name{getTipdata.Rd}
\alias{getTipdata}

\title{ Add species richness data to tips of phylogenetic tree  }
\description{
  Adds species richness (taxonomic) data to phylogenetic tree
}    
\usage{
getTipdata(tipdata, phy)

}


\arguments{
  \item{tipdata}{ either of the following: a vector of species diversities, where the name of each element corresponds to the name of a terminal in a phylogenetic tree, or a data frame with a single column, where the row names are the names of the terminals}
  \item{phy}{ a class 'phylo' phylogenetic tree, perhaps obtained by \code{read.tree}}

}
\details{
 \code{getTipdata} takes species richness data and adds a new component \code{phenotype} to your phylogenetic tree. The phenotype data could in principle be anything, but the function is designed to add species richness data to the tree to facilitate analysis by \code{fitLambda_1rate} and \code{fitNDR_2rate}. 
 
 This function requires that the names in \code{tipdata} exactly match \code{phy$tip.label}; if they don't, you will (probably) receive an error message. If you fail to supply a vector with names or a dataframe with row.names, the input vector will be matched exactly to the to tips in the order in which they are stored in your tree (e.g., \code{tipdata[1]} will correspond to \code{phy$tip.label[1]}). This is not recommended.
 
}
\value{
  a phylogenetic tree with the standard components of class 'phylo' plus the following component:
  \item{phenotype}{ Species diversities associated with each terminal taxon }
  
}
\references{ 

Rabosky, D. L., S. C. Donnellan, A. L. Talaba, and I. J. Lovette. 2007. Exceptional among-lineage variation in diversification rates during the radiation of Australia's largest vertebrate clade. Proc. Roy. Soc. Lond. Ser. B 274:2915-2923.


}
\author{ Dan Rabosky \email{DLR32@cornell.edu} }

\seealso{ \code{\link{fitNDR_2rate}}, \code{\link{getTipdata}}, \code{\link{lambda.stem.ml}}, \code{\link{fitNDR_1rate}}  }
\examples{

}
\keyword{ models }
