`yuleint2` <-
function(x, st1, st2)
{
  nvec <- 2:(length(x)+1)
  nv <- x[(x < st1) & (x >= st2)]
  lo <- max(nvec[x >= st1])
  up <- max(nvec[x >= st2])
  
  res <- list()
  if (st1 <= x[1])
    {nv <- c(st1, nv) - st2}
  else
    {nv <- nv - st2}

  smax <- (up-lo)/(lo*nv[1] + sum(nv[2:(up-lo +1)]))

  s1 <- sum(log(lo:(up-1)))
  s2 <- (up-lo)*log(smax)
  s3 <- -(lo*nv[1] + sum(nv[2:(up-lo+1)]))*smax
  res$smax <- smax
  res$LH <- s1 + s2 + s3

  return(res)

}

