`IDDX` <-
function(x)
{

  N <- length(x)+1
  b <- sort(x)
  z <- rev(c(b[1], diff(b)))
  #initial parameterization for nlm fx: estimates under yule model
  s1 <- IpureBirth(x)
  #print(s1$r1)
  res <- list()
  temp <- list()
  ddfunc <- function(r, v)
  {
    -(sum(log(2:(N-1))) + (N-2)*log(r) - sum(v*log(2:(N-1)))
    - sum(((2:N)^(1-v))*z *r) )
  }
  temp <- nlm(function(p) ddfunc(p[1], p[2]), c(s1$r1, 0), hessian = TRUE)

  res$LH <- -temp$minimum
  res$r1 <- temp$estimate[1]
  res$xp <- temp$estimate[2]

   return(res)
}

