`pureBirth` <-
function(x)
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

  return(res)
}

