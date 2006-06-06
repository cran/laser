getBtimes <- function(string = NULL, file = NULL)
{
  #requires 'ape' package.
  require(ape)
  if (is.null(string) && !is.null(file))
    tree <- read.tree(file)
  else if (!is.null(string) && is.null(file))
    tree <- read.tree(text = paste(string))
  else 
    stop("you must enter a filename or a character string\n")    
  if (!is.ultrametric(tree))
    stop("Tree is not ultrametric!")
  btimes <- branching.times(tree)
  btimes <- rev(sort(btimes))
  
  #strip node labels from btimes, as returned by 'branching.times'
  btimes <- as.numeric(btimes)
  return(btimes)

}
