`IpureBirth` <-
function(x)
{
  res <- list()
  temp <- yuleint2(x, x[1], 0)
  res$LH <- temp$LH
  res$r1 <- temp$smax
  return(res)
}

