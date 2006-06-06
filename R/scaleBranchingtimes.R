### scaleBranchingtimes
###
### Calibrates Branching Times to Basal Divergence Time
###

scaleBranchingtimes <- function(x, basal = 100)
{
  x <- as.numeric(x)
  x <- rev(sort(x))

  # 'basal' is the basal divergence time...
  scalefactor <- basal/x[1]
  scaledvec <- x*scalefactor
  return(scaledvec)
}
###########################