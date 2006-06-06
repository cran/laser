### pureBirth
###
### Fits Pure Birth Model to Set of Branching Times
###

pureBirth <- function(x)
{
  # calculates ML estimates of speciation rate
  #under yule model & associated LH of full set of branching times
  if (!is.numeric(x)) stop("object x not of class 'numeric'")
  res <- list()
  x <- rev(sort(x))
  temp <- yuleint2(x, x[1], 0)
  res$LH <- temp$LH
  res$aic <- (-2*res$LH) + 2
  res$r1 <- temp$smax
  cat("-----------------------------\n")
  cat("Model: pure birth\n")
  cat("Log-likelihood and parameter:\nLH", res$LH, "\nAIC", res$aic, "\nr", res$r1, "\n")

  return(res)
}