########################################################
### getBtimes.batch
###
### Get Branching Times For Batch of Phylogenies
###
###

getBtimes.batch <- function(file = NULL, format = "newick", basal = NULL)
{
  #if format = "newick", reads trees & converts.
  # if ="branchingtimes", simply read in btimes
  # see documentation for examples of correctly formatted files
  #needs package 'ape' to run correctly if format "newick"

  tset <- data.frame()
  if (format == "newick")
  {
    require(ape)
    if (!is.null(file))
      v <- read.tree(file) 
    else 
      stop("you must enter a filename where trees are stored\n")    
    
   
    Ntrees <- length(v)
    for (i in 1:Ntrees)
    {
      if ((i %% 10) == 0)
        cat("Processing tree", i, "\n")
      temp <- eval(parse(text = paste("branching.times(v$tree", i, ")", sep = "")))
      temp <- as.vector(temp)
      temp <- rev(sort(temp))
      #scales branching times by common basal divergence time if specified...#
      if (!is.null(basal))
      {
        scalefactor <- basal/temp[1]
        temp <- temp*scalefactor
      }
      if (i == 1)
      {
        tset <- temp
      }
      else
        tset <- rbind(tset, temp)
    }
  }
  else if (format == "branchingtimes")
  {
    tempvec <- scan(file)
    ntaxa <- tempvec[1]
    nsets <- tempvec[2]
    tempvec <- tempvec[3:length(tempvec)]
    tset <- matrix(tempvec, nrow = (ntaxa-1), ncol = nsets)
    tset <- t(tset)
  }
  return(tset)
}
#######################################################